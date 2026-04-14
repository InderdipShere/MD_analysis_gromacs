#!/usr/bin/env python3
"""
Calculates bond orientation with respect to a reference direction from GROMACS trajectories.

This script extends bond analysis to include orientation calculations by:
- Computing a reference direction vector for each molecule (from direction_species atoms)
- Detecting bonds between specified groups
- Computing the angle between each bond vector and the reference direction
- Filtering: only bonds where reference atom is from the same molecule as direction atoms
- Outputting average orientation angles and optional distribution histograms

Computational Details:
- Frame-by-frame trajectory processing with bond detection and angle calculation
- Periodic boundary condition handling via minimum image convention (MIC)
- Support for atom-based and center-of-mass (COM) reference/selection positions
- Vectorized distance and angle calculations using NumPy
- Bond detection based on cutoff distance (rcut)
- Orientation angle: arccos(bond_vec · direction_vec / (|bond_vec| * |direction_vec|))
- Multiple pair support with explicit pairing (ref[i] with sel[i])
"""
import argparse
import sys
import os
import logging
import time
import socket
import platform
import numpy as np


def get_cli_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-s', '--structure_file', required=True, help='Structure file (.gro) to initialize the system.')
    parser.add_argument('-f', '--traj_file', required=True, help='Trajectory file (.gro format).')
    parser.add_argument('-n', '--index_file', required=True, help='Index file (.ndx).')
    parser.add_argument('-o', '--output_file', default='orientation.xvg', help='Output file name for orientation statistics (default: orientation.xvg).')
    parser.add_argument('--ref', required=True, nargs='+', help='Name(s) of the reference group(s) in the index file.')
    parser.add_argument('--sel', required=True, nargs='+', help='Name(s) of the selection group(s) in the index file.')
    parser.add_argument('--ref_mol', type=str2bool, nargs='*', default=None, help='Use center-of-mass (COM) for reference group(s).')
    parser.add_argument('--sel_mol', type=str2bool, nargs='*', default=None, help='Use center-of-mass (COM) for selection group(s).')
    parser.add_argument('-direction_species', '--direction_species', required=True, nargs='+', help='Two atom names defining direction vector for each pair (same length as --ref).')
    parser.add_argument('--direction_mol', type=str2bool, nargs='*', default=None, help='Molecule-based direction (1) or atom-based (0). Default: False.')
    parser.add_argument('--mass_file', type=str, default=None, help='File containing atomic masses for COM calculation (required if using COM).')
    parser.add_argument('--rcut', type=float, nargs='+', required=True, help='Cutoff distance(s) (nm) for bond detection. Can be: single value (for all pairs) or multiple values (one per pair).')
    parser.add_argument('--begin', type=int, default=0, help='First frame to read from trajectory.')
    parser.add_argument('--end', type=int, default=-1, help='Last frame to read from trajectory. Use -1 for end.')
    parser.add_argument('--skip', type=int, default=1, help='Read every Nth frame.')
    parser.add_argument('-dist', '--dist', type=int, default=None, help='Generate orientation distribution histogram with N bins (0-180°). Output to output_dist.xvg')
    parser.add_argument('-debug', '--debug', action='store_true', help="Enable debug mode")
    
    return parser.parse_args()


def str2bool(v):
    """Convert string to boolean."""
    if isinstance(v, bool):
        return v
    v = v.lower()
    if v in ('1', 'true', 't', 'yes', 'y'):
        return True
    elif v in ('0', 'false', 'f', 'no', 'n'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected: 1/0, true/false, yes/no')


def parse_ndx(ndx_file):
    """Parse a GROMACS index file (.ndx) and return a dictionary of groups."""
    groups = {}
    current_group = None
    with open(ndx_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('[') and line.endswith(']'):
                current_group = line[1:-1].strip()
                groups[current_group] = []
            elif line and current_group:
                groups[current_group].extend(map(int, line.split()))
    # Convert to 0-based numpy arrays
    for name, indices in groups.items():
        groups[name] = np.array(indices, dtype=int) - 1
    return groups


def parse_mass_file(mass_file):
    """Parse mass file with atomic masses."""
    masses = {}
    current_mol = None
    with open(mass_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(';'):
                continue
            if line.startswith('[') and line.endswith(']'):
                current_mol = line[1:-1].strip()
                masses[current_mol] = {}
            elif '=' in line and current_mol is not None:
                atom, mass = line.split('=')
                masses[current_mol][atom.strip()] = float(mass)
    return masses


def read_gro_frame(f, skip=False):
    """
    Read a frame from a GROMACS .gro file stream.
    If skip=True, efficiently skip the frame without parsing.
    
    Returns:
        (coords, box, atom_info) or (None, None, None) if EOF or skip=True
    """
    title = f.readline()
    if not title:
        return None, None, None
    
    try:
        num_atoms_str = f.readline()
        if not num_atoms_str:
            return None, None, None
        num_atoms = int(num_atoms_str.strip())
        
        if skip:
            # Skip all atom lines and box line
            for _ in range(num_atoms + 1):
                f.readline()
            return None, None, None

        coords = np.zeros((num_atoms, 3), dtype=np.float32)
        atom_info = []
        for i in range(num_atoms):
            line = f.readline()
            # Extract: molno(0:5), resname(5:10), atomname(10:15), atomindex(15:20)
            atom_info.append((line[0:5].strip(), line[5:10].strip(), line[10:15].strip(), line[15:20].strip()))
            coords[i] = [float(line[20:28]), float(line[28:36]), float(line[36:44])]
        
        box_line = f.readline()
        box = np.array([float(x) for x in box_line.split()], dtype=np.float32)
        
        return coords, box, atom_info
    
    except (IOError, ValueError) as e:
        print(f"Error reading frame: {e}", file=sys.stderr)
        return None, None, None


def get_com(coords, groups_mol, mol_name, box=None):
    """
    Calculate center-of-mass for molecules in a group.
    
    Applies minimum image convention (MIC) to atomic coordinates within each molecule
    before calculating COM.
    
    Args:
        coords: (N_atoms, 3) coordinate array
        groups_mol: dict of molecular groups with mass info
        mol_name: name of the group to compute COM for
        box: (3,) box dimensions
    
    Returns:
        (N_mols, 3) array of COM coordinates
    """
    n_mols = len(groups_mol[mol_name])
    com = np.zeros((n_mols, 3), dtype=np.float32)
    
    for mol_id, mol in enumerate(groups_mol[mol_name]):
        indices = np.array(mol['atom_indices'], dtype=int)
        coords_mol = coords[indices].copy()
        masses = np.array(mol['atom_masses'], dtype=np.float32)
        
        # Apply MIC to atomic coordinates within molecule
        if box is not None and len(indices) > 1:
            ref_pos = coords_mol[0]
            for i in range(1, len(coords_mol)):
                delta = coords_mol[i] - ref_pos
                delta -= box[:3] * np.round(delta / box[:3])
                coords_mol[i] = ref_pos + delta
        
        # Calculate COM
        com[mol_id] = np.sum(coords_mol * masses[:, np.newaxis], axis=0) / mol['total_mass']
    
    return com


def compute_angle_degrees(vec1, vec2):
    """
    Compute angle in degrees between two vectors (or array of vectors).
    
    Args:
        vec1: (3,) or (N, 3) vector(s)
        vec2: (3,) or (N, 3) vector(s)
    
    Returns:
        angle(s) in degrees [0, 180]
    """
    # Normalize vectors
    norm1 = np.linalg.norm(vec1, axis=-1, keepdims=True)
    norm2 = np.linalg.norm(vec2, axis=-1, keepdims=True)
    
    # Avoid division by zero
    norm1[norm1 == 0] = 1.0
    norm2[norm2 == 0] = 1.0
    
    vec1_norm = vec1 / norm1
    vec2_norm = vec2 / norm2
    
    # Compute dot product
    if vec1.ndim == 1:
        cos_angle = np.dot(vec1_norm, vec2_norm)
    else:
        cos_angle = np.sum(vec1_norm * vec2_norm, axis=-1)
    
    # Clamp to [-1, 1] to handle numerical errors
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    
    # Convert to degrees
    angle_rad = np.arccos(cos_angle)
    angle_deg = np.degrees(angle_rad)
    
    return angle_deg


def compute_bonds_and_orientations_vectorized(ref_pos, sel_pos, direction_vec, box, rcut, 
                                             calculate_angles=True, return_angles=False, debug_first_ref=False):
    """
    Detect bonds and compute orientation angles between bond vectors and direction vector.
    
    Args:
        ref_pos: (N_ref, 3) array of reference positions
        sel_pos: (N_sel, 3) array of selection positions
        direction_vec: (3,) direction vector
        box: (3,) box dimensions
        rcut: cutoff distance for bond detection
        calculate_angles: if True, compute angles for all bonds
        return_angles: if True, return all bonded angles
        debug_first_ref: if True, return angle details for first ref particle
    
    Returns:
        Depending on flags:
        - n_bonds, avg_angle (if calculate_angles=True, not debug)
        - n_bonds, avg_angle, angle_array (if return_angles=True)
        - n_bonds, avg_angle, bonded_sel_indices, bonded_angles (if debug_first_ref=True)
    """
    # Vectorized distance calculation
    delta = sel_pos[np.newaxis, :, :] - ref_pos[:, np.newaxis, :]
    
    # Apply minimum image convention
    box_inv = 1.0 / box[:3]
    delta -= box[:3][np.newaxis, np.newaxis, :] * np.round(delta * box_inv[np.newaxis, np.newaxis, :])
    
    # Compute distances
    dist = np.linalg.norm(delta, axis=2)
    
    # Find bonds within cutoff
    bonds = dist <= rcut
    n_bonds = np.sum(bonds)
    
    if n_bonds == 0:
        if debug_first_ref:
            return n_bonds, 0.0, np.array([], dtype=int), np.array([], dtype=np.float32)
        elif return_angles:
            return n_bonds, 0.0, np.array([], dtype=np.float32)
        else:
            return n_bonds, 0.0
    
    if calculate_angles:
        # Compute angles for all bonded pairs
        bonded_delta = delta[bonds]  # (n_bonds, 3)
        
        # Compute angle with direction vector for each bond
        angles = compute_angle_degrees(bonded_delta, direction_vec)  # (n_bonds,)
        
        if debug_first_ref:
            # Get bonds for first ref particle
            first_ref_bonds = bonds[0, :]
            bonded_sel_indices = np.where(first_ref_bonds)[0]
            first_ref_angles = compute_angle_degrees(delta[0, bonded_sel_indices], direction_vec)
            avg_angle = np.mean(angles) if n_bonds > 0 else 0.0
            return n_bonds, avg_angle, bonded_sel_indices, first_ref_angles
        
        avg_angle = np.mean(angles)
        
        if return_angles:
            return n_bonds, avg_angle, angles
        else:
            return n_bonds, avg_angle
    
    else:
        # Just return bond count
        return n_bonds


def process_frame_positions(args, groups, groups_mol, coords, box, atom_info):
    """
    Extract and process positions for a single frame.
    Also extract direction vectors from direction_species.
    
    Direction vector for a molecule: vector from atom1 to atom2 in same molecule.
    atom2 can come from a different molecule (only needs to be within rcut).
    
    Returns:
        Dictionary with reference, selection positions and direction vectors for all pairs
    """
    frame_data = {}
    atom_info_array = np.array(atom_info)
    
    # Process reference groups
    for i, ref_name in enumerate(args.ref):
        ref_indices = groups[ref_name]
        if args.ref_mol[i]:
            ref_pos = get_com(coords, groups_mol, ref_name, box=box)
        else:
            ref_pos = coords[ref_indices]
        frame_data[('ref', i)] = ref_pos
    
    # Process selection groups
    for j, sel_name in enumerate(args.sel):
        sel_indices = groups[sel_name]
        if args.sel_mol[j]:
            sel_pos = get_com(coords, groups_mol, sel_name, box=box)
        else:
            sel_pos = coords[sel_indices]
        frame_data[('sel', j)] = sel_pos
    
    # Process direction vectors for each pair
    for i, pair_direction in enumerate(args.direction_species_pairs):
        atom1_name, atom2_name = pair_direction
        
        # Get all indices for these atom names
        atom1_indices = np.where(atom_info_array[:, 2] == atom1_name)[0]
        atom2_indices = np.where(atom_info_array[:, 2] == atom2_name)[0]
        
        if len(atom1_indices) == 0 or len(atom2_indices) == 0:
            logging.warning(f"Pair {i}: Could not find atoms {atom1_name} or {atom2_name}")
            frame_data[('dir_vec', i)] = None
            continue
        
        # For each molecule, compute direction vector from atom1 in that molecule to atom2 (any molecule)
        direction_vectors = {}
        
        # For each atom1 occurrence (each molecule that has atom1)
        for idx1 in atom1_indices:
            mol_no = atom_info[idx1][0]  # Molecule ID where atom1 is located
            
            # Find closest atom2 (could be from any molecule)
            pos1 = coords[idx1]
            atom2_positions = coords[atom2_indices]
            
            # Apply MIC for distance calculation
            deltas = atom2_positions - pos1
            deltas -= box[:3] * np.round(deltas / box[:3])
            distances = np.linalg.norm(deltas, axis=1)
            
            # Find closest atom2
            closest_idx2 = atom2_indices[np.argmin(distances)]
            pos2 = coords[closest_idx2]
            
            # Compute direction vector with MIC
            delta = pos2 - pos1
            delta -= box[:3] * np.round(delta / box[:3])
            
            direction_vectors[mol_no] = delta
        
        frame_data[('dir_vec', i)] = direction_vectors
    
    return frame_data


def calculate_orientations(args, groups, box, groups_mol, atom_info):
    """
    Calculate bond orientation statistics from trajectory.
    
    For each reference atom, determine its molecule ID and use that molecule's
    direction vector to compute bond angles.
    
    Returns:
        (orientations_per_frame, total_frames, orientation_angles)
        orientations_per_frame: (n_pairs, n_frames) array of average angles
        orientation_angles: list of angle arrays if --dist specified, else None
    """
    npairs = len(args.ref)
    frame_orientations = []
    orientation_angles = [[] for _ in range(npairs)] if args.dist else None
    total_frames = 0
    
    with open(args.traj_file, 'r') as f:
        # Skip to begin frame
        for _ in range(args.begin):
            read_gro_frame(f, skip=True)
        
        if args.end == -1:
            # Process until EOF
            while True:
                coords, box_frame, atom_info_frame = read_gro_frame(f, skip=False)
                if coords is None:
                    break
                
                # Process frame
                frame_data = process_frame_positions(args, groups, groups_mol, coords, box_frame, atom_info_frame)
                
                frame_orientations_data = []
                
                for i, (ref_name, sel_name) in enumerate(zip(args.ref, args.sel)):
                    ref_indices = groups[ref_name]
                    sel_pos = frame_data[('sel', i)]
                    direction_vectors = frame_data[('dir_vec', i)]
                    rcut_pair = args.rcut[i]
                    
                    if direction_vectors is None or len(direction_vectors) == 0:
                        frame_orientations_data.append(0.0)
                        continue
                    
                    # Process each reference atom individually
                    all_angles_for_frame = []
                    
                    for ref_idx in ref_indices:
                        mol_no = atom_info_frame[ref_idx][0]
                        
                        # Check if this molecule has a direction vector
                        if mol_no not in direction_vectors:
                            continue
                        
                        direction_vec = direction_vectors[mol_no]
                        ref_pos = coords[ref_idx:ref_idx+1]  # Shape (1, 3)
                        
                        # Compute bonds and angles for this reference atom
                        if args.dist:
                            n_bonds, avg_angle, angles = compute_bonds_and_orientations_vectorized(
                                ref_pos, sel_pos, direction_vec, box_frame, rcut_pair,
                                calculate_angles=True, return_angles=True
                            )
                            if len(angles) > 0:
                                all_angles_for_frame.extend(angles)
                                orientation_angles[i].extend(angles)
                        else:
                            n_bonds, avg_angle = compute_bonds_and_orientations_vectorized(
                                ref_pos, sel_pos, direction_vec, box_frame, rcut_pair,
                                calculate_angles=True, return_angles=False
                            )
                            if n_bonds > 0:
                                all_angles_for_frame.append(avg_angle)
                    
                    # Average angles across all processed reference atoms in this frame
                    if len(all_angles_for_frame) > 0:
                        frame_avg_angle = np.mean(all_angles_for_frame)
                    else:
                        frame_avg_angle = 0.0
                    
                    frame_orientations_data.append(frame_avg_angle)
                
                frame_orientations.append(frame_orientations_data)
                total_frames += 1
                
                # Skip frames based on skip parameter
                for _ in range(args.skip - 1):
                    read_gro_frame(f, skip=True)
        else:
            # Process specific frame range
            frame_idx = args.begin
            while frame_idx < args.end:
                coords, box_frame, atom_info_frame = read_gro_frame(f, skip=False)
                if coords is None:
                    logging.warning(f"Reached EOF before frame {args.end}")
                    break
                
                frame_data = process_frame_positions(args, groups, groups_mol, coords, box_frame, atom_info_frame)
                frame_orientations_data = []
                
                for i, (ref_name, sel_name) in enumerate(zip(args.ref, args.sel)):
                    ref_indices = groups[ref_name]
                    sel_pos = frame_data[('sel', i)]
                    direction_vectors = frame_data[('dir_vec', i)]
                    rcut_pair = args.rcut[i]
                    
                    if direction_vectors is None or len(direction_vectors) == 0:
                        frame_orientations_data.append(0.0)
                        continue
                    
                    all_angles_for_frame = []
                    
                    for ref_idx in ref_indices:
                        mol_no = atom_info_frame[ref_idx][0]
                        
                        if mol_no not in direction_vectors:
                            continue
                        
                        direction_vec = direction_vectors[mol_no]
                        ref_pos = coords[ref_idx:ref_idx+1]
                        
                        if args.dist:
                            n_bonds, avg_angle, angles = compute_bonds_and_orientations_vectorized(
                                ref_pos, sel_pos, direction_vec, box_frame, rcut_pair,
                                calculate_angles=True, return_angles=True
                            )
                            if len(angles) > 0:
                                all_angles_for_frame.extend(angles)
                                orientation_angles[i].extend(angles)
                        else:
                            n_bonds, avg_angle = compute_bonds_and_orientations_vectorized(
                                ref_pos, sel_pos, direction_vec, box_frame, rcut_pair,
                                calculate_angles=True, return_angles=False
                            )
                            if n_bonds > 0:
                                all_angles_for_frame.append(avg_angle)
                    
                    if len(all_angles_for_frame) > 0:
                        frame_avg_angle = np.mean(all_angles_for_frame)
                    else:
                        frame_avg_angle = 0.0
                    
                    frame_orientations_data.append(frame_avg_angle)
                
                frame_orientations.append(frame_orientations_data)
                total_frames += 1
                frame_idx += args.skip
                
                for _ in range(args.skip - 1):
                    read_gro_frame(f, skip=True)
    
    # Convert to numpy arrays
    orientations_per_frame = np.array(frame_orientations).T  # Shape: (n_pairs, n_frames)
    
    # Convert orientation_angles lists to numpy arrays
    if args.dist:
        orientation_angles = [np.array(d, dtype=np.float32) if d else None for d in orientation_angles]
    
    return orientations_per_frame, total_frames, orientation_angles


def normalize_and_output(args, orientations_per_frame, total_frames, ref_counts, sel_counts):
    """
    Normalize and write orientation statistics to file(s).
    
    Outputs average orientation (angle in degrees) per frame for each pair.
    """
    if total_frames == 0:
        logging.warning("No frames were processed. Exiting.")
        return
    
    # Create frame numbers for output
    frame_numbers = np.arange(total_frames)
    
    # Orientation output
    output_data = [frame_numbers]
    header = "# frame"
    
    logging.info("--- Orientation Output ---")
    logging.info(f"Total frames processed: {total_frames}")
    
    # Log rcut and direction species
    logging.info("\n--- Configuration ---")
    if len(set(args.rcut)) == 1:
        logging.info(f"Cutoff distance (rcut): {args.rcut[0]:.3f} nm (all pairs)")
    else:
        logging.info(f"Cutoff distances (rcut) per pair:")
        for i, (ref_name, sel_name) in enumerate(zip(args.ref, args.sel)):
            logging.info(f"  {ref_name}-{sel_name}: {args.rcut[i]:.3f} nm")
    
    logging.info("\nDirection vectors:")
    for i, direction_pair in enumerate(args.direction_species_pairs):
        logging.info(f"  Pair {i}: {direction_pair[0]} -> {direction_pair[1]}")
    
    for i, (ref_name, sel_name) in enumerate(zip(args.ref, args.sel)):
        N_ref = int(ref_counts[i])
        N_sel = int(sel_counts[i])
        
        logging.info(f"\nOrientation statistics for pair: {ref_name}-{sel_name} (reference direction: {args.direction_species_pairs[i][0]}-{args.direction_species_pairs[i][1]})")
        logging.info(f"  N_ref ({ref_name}): {N_ref}")
        logging.info(f"  N_sel ({sel_name}): {N_sel}")
        
        # Orientation angles per frame
        angles = orientations_per_frame[i, :]
        valid_angles = angles[angles > 0]  # Filter out frames with no bonds
        
        if len(valid_angles) > 0:
            logging.info(f"  Mean orientation angle: {np.mean(valid_angles):.2f}°")
            logging.info(f"  Min orientation angle: {np.min(valid_angles):.2f}°")
            logging.info(f"  Max orientation angle: {np.max(valid_angles):.2f}°")
            logging.info(f"  Std dev: {np.std(valid_angles):.2f}°")
        else:
            logging.warning(f"  No bonds found for {ref_name}-{sel_name}")
        
        output_data.append(angles)
        header += f" Orientation_{ref_name}-{sel_name}(deg)"
    
    # Write orientation output
    with open(args.output_file, 'w') as f:
        f.write(header + "\n")
        np.savetxt(f, np.column_stack(output_data), fmt='%10.6f')
        f.write("# Orientation: angle between bond vector (ref->sel) and direction vector\n")
        f.write("# Total number of Ref:\n")
        f.write(f"# {', '.join([f'{name}: {int(count)}' for name, count in zip(args.ref, ref_counts)])}\n")
        f.write("# Total number of sel:\n")
        f.write(f"# {', '.join([f'{name}: {int(count)}' for name, count in zip(args.sel, sel_counts)])}\n")
    
    logging.info(f"\nOrientation data written to {args.output_file}")


def generate_orientation_distribution(args, orientation_angles, ref_counts, sel_counts):
    """
    Generate and write orientation distribution as normalized histograms.
    
    The histogram is normalized such that:
    Integral of histogram over angle range = total number of bonds considered
    
    Args:
        args: command-line arguments (contains -dist N)
        orientation_angles: (n_pairs,) list of angle arrays
        ref_counts: reference particle/molecule counts
        sel_counts: selection particle/molecule counts
    """
    if args.dist is None or orientation_angles is None:
        return
    
    logging.info("\n--- Orientation Distribution Histogram (0-180°) ---")
    
    n_pairs = len(args.ref)
    n_bins = args.dist
    
    # Create output file
    dist_output_file = args.output_file.split('.')[0] + '_dist.' + args.output_file.split('.')[-1]
    
    # Create bin edges: from 0 to 180 degrees
    bin_edges = np.linspace(0, 180, n_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    bin_width = bin_edges[1] - bin_edges[0]
    
    output_data = [bin_centers]
    header = "# angle(degrees)"
    
    for i, (ref_name, sel_name) in enumerate(zip(args.ref, args.sel)):
        if orientation_angles[i] is None or len(orientation_angles[i]) == 0:
            count_dist = np.zeros(n_bins)
            logging.info(f"Orientation distribution for pair {ref_name}-{sel_name}: No angles found")
        else:
            # Create histogram
            hist, _ = np.histogram(orientation_angles[i], bins=bin_edges)
            
            # Normalize histogram so that integral = total number of bonds
            # This gives: sum(count_dist) * bin_width = total_count
            total_count = np.sum(hist)
            count_dist = hist / bin_width if bin_width > 0 else hist
            
            logging.info(f"Orientation distribution for pair {ref_name}-{sel_name}:")
            logging.info(f"  Total bonded pairs: {len(orientation_angles[i])}")
            logging.info(f"  Min angle: {np.min(orientation_angles[i]):.2f}°")
            logging.info(f"  Max angle: {np.max(orientation_angles[i]):.2f}°")
            logging.info(f"  Mean angle: {np.mean(orientation_angles[i]):.2f}°")
            logging.info(f"  Median angle: {np.median(orientation_angles[i]):.2f}°")
            logging.info(f"  Std dev: {np.std(orientation_angles[i]):.2f}°")
            logging.info(f"  Histogram bins: {n_bins}, range: 0-180°, bin width: {bin_width:.2f}°")
            logging.info(f"  Integration check (should equal {total_count}): {np.sum(count_dist) * bin_width:.2f}")
        
        output_data.append(count_dist)
        header += f" Hist_{ref_name}-{sel_name}"
    
    # Write distribution output
    with open(dist_output_file, 'w') as f:
        f.write(header + "\n")
        np.savetxt(f, np.column_stack(output_data), fmt='%12.6f')
        f.write("# Normalized histogram: H(θ) = count(θ bin) / bin_width\n")
        f.write("# Integral of H(θ) * dθ over angles 0-180° = total number of bonds\n")
        f.write("# Angle θ is between bond vector (ref->sel) and reference direction vector\n")
        f.write("# Total number of Ref:\n")
        f.write(f"# {', '.join([f'{name}: {int(count)}' for name, count in zip(args.ref, ref_counts)])}\n")
        f.write("# Total number of sel:\n")
        f.write(f"# {', '.join([f'{name}: {int(count)}' for name, count in zip(args.sel, sel_counts)])}\n")
    
    logging.info(f"Orientation distribution written to {dist_output_file}")


def main():
    """Main function to run the bond orientation analysis from GROMACS trajectory."""
    t_start = time.time()
    
    args = get_cli_args()
    need_mass_file = False

    # Validate pair structure
    nref = len(args.ref)
    nsel = len(args.sel)
    ndir = len(args.direction_species)
    
    if nref != nsel:
        raise ValueError(
            f"--ref and --sel must have the same length. "
            f"Got {nref} ref groups and {nsel} sel groups.")
    
    if ndir != nref * 2:
        raise ValueError(
            f"--direction_species must have 2*{nref} = {nref*2} elements (2 atoms per pair). "
            f"Got {ndir} elements.")
    
    # Parse direction species pairs
    args.direction_species_pairs = [(args.direction_species[2*i], args.direction_species[2*i+1]) for i in range(nref)]
    
    npairs = nref
    
    # Validate and process rcut values
    if len(args.rcut) == 1:
        args.rcut = args.rcut * npairs
    elif len(args.rcut) != npairs:
        raise ValueError(
            f"--rcut must have length 1 or {npairs}, but got {len(args.rcut)}")
    
    # Process ref_mol flags
    if args.ref_mol is None:
        args.ref_mol = [False] * nref
    elif len(args.ref_mol) == 1:
        if args.ref_mol[0]:
            need_mass_file = True
        args.ref_mol = args.ref_mol * nref
    elif len(args.ref_mol) != nref:
        raise ValueError(f"--ref_mol must have length 1 or {nref}, but got {len(args.ref_mol)}")
    else:
        need_mass_file = True

    # Process sel_mol flags
    if args.sel_mol is None:
        args.sel_mol = [False] * nsel
    elif len(args.sel_mol) == 1:
        if args.sel_mol[0]:
            need_mass_file = True
        args.sel_mol = args.sel_mol * nsel
    elif len(args.sel_mol) != nsel:
        raise ValueError(f"--sel_mol must have length 1 or {nsel}, but got {len(args.sel_mol)}")
    else:
        need_mass_file = True
    
    # Process direction_mol flags
    if args.direction_mol is None:
        args.direction_mol = [False] * nref
    elif len(args.direction_mol) == 1:
        if args.direction_mol[0]:
            need_mass_file = True
        args.direction_mol = args.direction_mol * nref
    elif len(args.direction_mol) != nref:
        raise ValueError(f"--direction_mol must have length 1 or {nref}, but got {len(args.direction_mol)}")
    else:
        need_mass_file = True
    
    if need_mass_file:
        if args.mass_file is None:
            raise ValueError("Mass file is required for COM calculations but not provided.")

    # Setup Logging
    log_file = os.path.splitext(args.output_file)[0] + '.log'
    logging.basicConfig(
        level=logging.INFO,
        format='%(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler(sys.stdout)
        ]
    )

    # Get system information
    hostname = socket.gethostname()
    platform_info = f"{platform.system()} {platform.release()}"
    python_version = platform.python_version()

    logging.info("=" * 70)
    logging.info("BOND ORIENTATION ANALYSIS FROM TRAJECTORY")
    logging.info("=" * 70)
    logging.info(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(f"Hostname: {hostname}")
    logging.info(f"Platform: {platform_info}")
    logging.info(f"Python: {python_version}")
    logging.info("")
    logging.info("Input Parameters:")
    for k, v in vars(args).items():
        if k != 'direction_species_pairs':  # Don't print derived field
            logging.info(f"  {k}: {v}")
    
    logging.info("\n--- Initialization ---")
    logging.info("Parsing index file...")
    groups = parse_ndx(args.index_file)
    
    if args.debug:
        for group_name in args.ref:
            atom_indices = groups[group_name] + 1
            logging.info(f"Reference group '{group_name}' has {len(groups[group_name])} atoms.")
        for group_name in args.sel:
            atom_indices = groups[group_name] + 1
            logging.info(f"Selection group '{group_name}' has {len(groups[group_name])} atoms.")
    
    if need_mass_file:
        logging.info("Getting masses for COM calculation...")
        masses = parse_mass_file(args.mass_file)
        if masses is None:
            raise ValueError("Mass file not parsed correctly.")

    logging.info("Initializing system from structure file...")    
    groups_mol = None
    with open(args.structure_file, 'r') as f:
        coords, box, atom_info = read_gro_frame(f)
        atom_info = np.array(atom_info)
        
        if need_mass_file:
            groups_mol = {}
            for group_name, atom_indices in groups.items():
                if group_name not in args.ref + args.sel:
                    continue
                
                mols_in_group = {}
                for idx in atom_indices:
                    mol_no = atom_info[idx][0]
                    atom_name = atom_info[idx][2]
                    
                    if mol_no not in mols_in_group:
                        mols_in_group[mol_no] = {
                            "mol_index": mol_no,
                            "atom_indices": [],
                            "atom_names": [],
                            "atom_masses": [],
                            "total_mass": 0.0
                        }
                    
                    mols_in_group[mol_no]["atom_indices"].append(idx)
                    mols_in_group[mol_no]["atom_names"].append(atom_name)
                    
                    try:
                        mass = masses[group_name][atom_name]
                    except KeyError:
                        mass = 0.0
                    
                    mols_in_group[mol_no]["atom_masses"].append(mass)
                    mols_in_group[mol_no]["total_mass"] = np.sum(mols_in_group[mol_no]["atom_masses"])
                
                groups_mol[group_name] = list(mols_in_group.values())
            
            logging.info(f"Groups prepared for COM calculation: {list(groups_mol.keys())}")

    if coords is None:
        logging.error("Issue with the structure file.")
        sys.exit(1)

    logging.info(f"Box dimensions: {box}")
    if len(set(args.rcut)) == 1:
        logging.info(f"Cutoff distance (rcut): {args.rcut[0]:.3f} nm (all pairs)")
    else:
        logging.info(f"Cutoff distances (rcut):")
        for i, (ref_name, sel_name) in enumerate(zip(args.ref, args.sel)):
            logging.info(f"  {ref_name}-{sel_name}: {args.rcut[i]:.3f} nm")
    
    # Get initial counts
    ref_counts = np.zeros(len(args.ref))
    sel_counts = np.zeros(len(args.sel))
    
    for i, ref_name in enumerate(args.ref):
        ref_indices = groups[ref_name]
        if args.ref_mol[i]:
            ref_counts[i] = len(groups_mol[ref_name])
            logging.info(f"No of molecules in ref group '{ref_name}': {ref_counts[i]:.0f}")
        else:
            ref_counts[i] = len(ref_indices)
            logging.info(f"No of atoms in ref group '{ref_name}': {ref_counts[i]:.0f}")
    
    for i, sel_name in enumerate(args.sel):
        sel_indices = groups[sel_name]
        if args.sel_mol[i]:
            sel_counts[i] = len(groups_mol[sel_name])
            logging.info(f"No of molecules in sel group '{sel_name}': {sel_counts[i]:.0f}")
        else:
            sel_counts[i] = len(sel_indices)
            logging.info(f"No of atoms in sel group '{sel_name}': {sel_counts[i]:.0f}")
    
    logging.info("Reference group counts: " + ", ".join([f"{name}: {int(count)}" for name, count in zip(args.ref, ref_counts)]))
    logging.info("Selection group counts: " + ", ".join([f"{name}: {int(count)}" for name, count in zip(args.sel, sel_counts)]))

    logging.info("\n--- Calculation ---")
    logging.info("Computing bond orientations using vectorized calculations...")
    
    orientations_per_frame, total_frames, orientation_angles = calculate_orientations(args, groups, box, groups_mol, atom_info)
    logging.info(f"Total frames processed: {total_frames}")
    
    logging.info("\n--- Output ---")
    normalize_and_output(args, orientations_per_frame, total_frames, ref_counts, sel_counts)
    
    # Generate orientation distribution if requested
    if args.dist:
        generate_orientation_distribution(args, orientation_angles, ref_counts, sel_counts)
    
    # Report timing information
    t_end = time.time()
    t_elapsed = t_end - t_start
    t_per_frame = t_elapsed / max(1, total_frames) if total_frames > 0 else 0.0
    
    logging.info("\n" + "=" * 70)
    logging.info("EXECUTION SUMMARY")
    logging.info("=" * 70)
    logging.info(f"Total execution time: {t_elapsed:.2f} seconds")
    if total_frames > 0:
        logging.info(f"Time per frame: {t_per_frame:.4f} seconds")
        logging.info(f"Frames per second: {total_frames / t_elapsed:.2f}")
    logging.info(f"Completed: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info("Bond orientation analysis complete!")
    logging.info("=" * 70)


if __name__ == '__main__':
    main()
