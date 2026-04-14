#!/usr/bin/env python3
"""
Calculates bond statistics from a GROMACS trajectory.

This script manually parses GROMACS .gro and .ndx files to detect bonds
between specified groups of atoms or molecules (using center-of-mass) and
calculate average number of bonds and optionally bond lengths.

Computational Details:
- Frame-by-frame trajectory processing with bond detection
- Periodic boundary condition handling via minimum image convention (MIC)
- Support for atom-based and center-of-mass (COM) reference/selection positions
- Vectorized distance calculations using NumPy broadcasting
- Bond detection based on cutoff distance (rcut)
- Optional bond length statistics
- Multi-pair output with support for multiple reference and selection groups
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
    parser.add_argument('-o', '--output_file', default='bonds.xvg', help='Output file name for bond statistics (default: bonds.xvg).')
    parser.add_argument('--BL', action='store_true', help='Also calculate average bond length.')
    parser.add_argument('--ref', required=True, nargs='+', help='Name(s) of the reference group(s) in the index file.')
    parser.add_argument('--sel', required=True, nargs='+', help='Name(s) of the selection group(s) in the index file.')
    parser.add_argument('--ref_mol', type=str2bool, nargs='*', default=None, help='Use center-of-mass (COM) for reference group(s).')
    parser.add_argument('--sel_mol', type=str2bool, nargs='*', default=None, help='Use center-of-mass (COM) for selection group(s).')
    parser.add_argument('--mass_file', type=str, default=None, help='File containing atomic masses for COM calculation (required if using COM).')
    parser.add_argument('--rcut', type=float, nargs='+', required=True, help='Cutoff distance(s) (nm) for bond detection. Can be: single value (for all pairs) or multiple values (one per pair).')
    parser.add_argument('--begin', type=int, default=0, help='First frame to read from trajectory.')
    parser.add_argument('--end', type=int, default=-1, help='Last frame to read from trajectory. Use -1 for end.')
    parser.add_argument('--skip', type=int, default=1, help='Read every Nth frame.')
    parser.add_argument('-dist', '--dist', type=int, default=None, help='Generate bond distance distribution histogram with N bins. Output to output_dist.xvg')
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
    before calculating COM. This ensures molecules split across box boundaries are
    handled correctly.
    
    Args:
        coords: (N_atoms, 3) coordinate array
        groups_mol: dict of molecular groups with mass info
        mol_name: name of the group to compute COM for
        box: (3,) box dimensions. If provided, MIC is applied to atomic coords within molecules.
    
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
        # Use first atom as reference to prevent molecule splitting across boundary
        if box is not None and len(indices) > 1:
            ref_pos = coords_mol[0]
            for i in range(1, len(coords_mol)):
                delta = coords_mol[i] - ref_pos
                delta -= box[:3] * np.round(delta / box[:3])
                coords_mol[i] = ref_pos + delta
        
        # Calculate COM
        com[mol_id] = np.sum(coords_mol * masses[:, np.newaxis], axis=0) / mol['total_mass']
    
    return com


def compute_bonds_vectorized(ref_pos, sel_pos, box, rcut, calculate_lengths=False, debug_first_ref=False, return_distances=False):
    """
    Detect bonds between reference and selection positions using vectorized operations.
    
    Uses NumPy broadcasting to compute all pairwise distances efficiently with
    minimum image convention (MIC). A bond is detected if distance <= rcut.
    
    Args:
        ref_pos: (N_ref, 3) array of reference positions
        sel_pos: (N_sel, 3) array of selection positions
        box: (3,) box dimensions
        rcut: cutoff distance for bond detection
        calculate_lengths: if True, also return bond lengths
        debug_first_ref: if True, return bond details for first ref particle
        return_distances: if True, also return all bonded distances
    
    Returns:
        If debug_first_ref=False:
            If calculate_lengths=False and return_distances=False:
                n_bonds: number of bonds detected
            If calculate_lengths=True and return_distances=False:
                (n_bonds, avg_bond_length): tuple with count and average length
            If calculate_lengths=False and return_distances=True:
                (n_bonds, bonded_distances): all distances of bonded pairs
            If both True:
                (n_bonds, avg_bond_length, bonded_distances)
        If debug_first_ref=True:
            (n_bonds, avg_bond_length, bond_indices, bond_distances) where:
                bond_indices: indices of sel atoms bonded with first ref
                bond_distances: distances of those bonds
    """
    # Vectorized distance calculation using broadcasting
    # delta shape: (N_ref, N_sel, 3)
    delta = sel_pos[np.newaxis, :, :] - ref_pos[:, np.newaxis, :]
    
    # Apply minimum image convention (MIC) - vectorized
    # Ensures we use the shortest distance considering periodic boundaries
    box_inv = 1.0 / box[:3]
    delta -= box[:3][np.newaxis, np.newaxis, :] * np.round(delta * box_inv[np.newaxis, np.newaxis, :])
    
    # Compute distances: shape (N_ref, N_sel)
    dist = np.linalg.norm(delta, axis=2)
    
    # Find bonds within cutoff (distance <= rcut)
    bonds = dist <= rcut
    n_bonds = np.sum(bonds)
    
    if debug_first_ref:
        # Get bonds for first ref particle
        first_ref_bonds = bonds[0, :]
        first_ref_distances = dist[0, :]
        bonded_sel_indices = np.where(first_ref_bonds)[0]
        bonded_distances = first_ref_distances[bonded_sel_indices]
        
        if calculate_lengths and n_bonds > 0:
            bond_lengths = dist[bonds]
            avg_bond_length = np.mean(bond_lengths)
            return n_bonds, avg_bond_length, bonded_sel_indices, bonded_distances
        elif calculate_lengths:
            return n_bonds, 0.0, bonded_sel_indices, bonded_distances
        else:
            return n_bonds, bonded_sel_indices, bonded_distances
    
    if return_distances:
        if n_bonds > 0:
            bonded_distances = dist[bonds]
        else:
            bonded_distances = np.array([], dtype=np.float32)
        
        if calculate_lengths and n_bonds > 0:
            avg_bond_length = np.mean(bonded_distances)
            return n_bonds, avg_bond_length, bonded_distances
        elif calculate_lengths:
            return n_bonds, 0.0, bonded_distances
        else:
            return n_bonds, bonded_distances
    
    if calculate_lengths and n_bonds > 0:
        bond_lengths = dist[bonds]
        avg_bond_length = np.mean(bond_lengths)
        return n_bonds, avg_bond_length
    elif calculate_lengths:
        return n_bonds, 0.0
    else:
        return n_bonds


def process_frame_positions(args, groups, groups_mol, coords, box):
    """
    Extract and process positions for a single frame.
    
    Applies MIC to molecular coordinates before calculating COM to ensure
    molecules split across box boundaries are handled correctly.
    
    Returns:
        Dictionary with reference and selection positions for all pairs
    """
    frame_data = {}
    
    for i, ref_name in enumerate(args.ref):
        ref_indices = groups[ref_name]
        if args.ref_mol[i]:
            ref_pos = get_com(coords, groups_mol, ref_name, box=box)
        else:
            ref_pos = coords[ref_indices]
        frame_data[('ref', i)] = ref_pos
    
    for j, sel_name in enumerate(args.sel):
        sel_indices = groups[sel_name]
        if args.sel_mol[j]:
            sel_pos = get_com(coords, groups_mol, sel_name, box=box)
        else:
            sel_pos = coords[sel_indices]
        frame_data[('sel', j)] = sel_pos
    
    return frame_data


def calculate_bonds(args, groups, box, groups_mol):
    """
    Calculate bond statistics from trajectory.
    
    Args:
        args: command-line arguments
        groups: dictionary of atom groups
        box: box dimensions from initial structure
        groups_mol: dictionary of molecular groups
    
    Returns:
        (bond_counts, bond_lengths_avg, total_frames, bond_distances)
        bond_counts: (n_pairs, n_frames) array
        bond_lengths_avg: (n_pairs, n_frames) array if --BL flag set, else None
        bond_distances: list of distance arrays if --dist specified, else None
    """
    npairs = len(args.ref)  # Explicit pairing: ref[i] with sel[i]
    frame_bonds = []
    frame_bond_lengths = [] if args.BL else None
    bond_distances = [[] for _ in range(npairs)] if args.dist else None
    total_frames = 0
    
    with open(args.traj_file, 'r') as f:
        # Skip to begin frame
        for _ in range(args.begin):
            read_gro_frame(f, skip=True)
        
        if args.end == -1:
            # Process until EOF
            while True:
                coords, box_frame, _ = read_gro_frame(f, skip=False)
                if coords is None:
                    break
                
                # Process frame positions with MIC for molecules
                frame_data = process_frame_positions(args, groups, groups_mol, coords, box_frame)
                
                # Compute bonds for all pairs (explicit pairing: ref[i] with sel[i])
                frame_bonds_data = []
                frame_lengths_data = [] if args.BL else None
                
                for i, (ref_name, sel_name) in enumerate(zip(args.ref, args.sel)):
                    ref_pos = frame_data[('ref', i)]
                    sel_pos = frame_data[('sel', i)]
                    rcut_pair = args.rcut[i]
                    
                    # Debug: print bonds for first ref molecule/atom
                    debug_flag = args.debug and i == 0 and total_frames == 0
                    
                    if args.BL and args.dist:
                        if debug_flag:
                            result = compute_bonds_vectorized(
                                ref_pos, sel_pos, box_frame, rcut_pair, calculate_lengths=True, debug_first_ref=True, return_distances=True
                            )
                            n_bonds, avg_length, bonded_sel_indices, bonded_distances = result[:4]
                            logging.info(f"\n[DEBUG] Frame {total_frames}, Pair {ref_name}-{sel_name}:")
                            logging.info(f"  First {ref_name}: {bonded_sel_indices.size} bonds")
                            for sel_idx, dist_val in zip(bonded_sel_indices, bonded_distances):
                                logging.info(f"    -> {sel_name}[{sel_idx}]: {dist_val:.6f} nm")
                        else:
                            n_bonds, avg_length, pair_distances = compute_bonds_vectorized(
                                ref_pos, sel_pos, box_frame, rcut_pair, calculate_lengths=True, return_distances=True
                            )
                            bond_distances[i].extend(pair_distances)
                        frame_bonds_data.append(int(n_bonds))
                        frame_lengths_data.append(float(avg_length))
                    elif args.BL:
                        if debug_flag:
                            result = compute_bonds_vectorized(
                                ref_pos, sel_pos, box_frame, rcut_pair, calculate_lengths=True, debug_first_ref=True
                            )
                            n_bonds, avg_length, bonded_sel_indices, bonded_distances = result
                            logging.info(f"\n[DEBUG] Frame {total_frames}, Pair {ref_name}-{sel_name}:")
                            logging.info(f"  First {ref_name}: {bonded_sel_indices.size} bonds")
                            for sel_idx, dist_val in zip(bonded_sel_indices, bonded_distances):
                                logging.info(f"    -> {sel_name}[{sel_idx}]: {dist_val:.6f} nm")
                        else:
                            n_bonds, avg_length = compute_bonds_vectorized(
                                ref_pos, sel_pos, box_frame, rcut_pair, calculate_lengths=True
                            )
                        frame_bonds_data.append(int(n_bonds))
                        frame_lengths_data.append(float(avg_length))
                    elif args.dist:
                        if debug_flag:
                            result = compute_bonds_vectorized(
                                ref_pos, sel_pos, box_frame, rcut_pair, calculate_lengths=False, debug_first_ref=True, return_distances=True
                            )
                            n_bonds, bonded_sel_indices, bonded_distances = result[:3]
                            logging.info(f"\n[DEBUG] Frame {total_frames}, Pair {ref_name}-{sel_name}:")
                            logging.info(f"  First {ref_name}: {bonded_sel_indices.size} bonds")
                            for sel_idx, dist_val in zip(bonded_sel_indices, bonded_distances):
                                logging.info(f"    -> {sel_name}[{sel_idx}]: {dist_val:.6f} nm")
                        else:
                            n_bonds, pair_distances = compute_bonds_vectorized(
                                ref_pos, sel_pos, box_frame, rcut_pair, calculate_lengths=False, return_distances=True
                            )
                            bond_distances[i].extend(pair_distances)
                        frame_bonds_data.append(int(n_bonds))
                    else:
                        if debug_flag:
                            result = compute_bonds_vectorized(
                                ref_pos, sel_pos, box_frame, rcut_pair, calculate_lengths=False, debug_first_ref=True
                            )
                            n_bonds, bonded_sel_indices, bonded_distances = result
                            logging.info(f"\n[DEBUG] Frame {total_frames}, Pair {ref_name}-{sel_name}:")
                            logging.info(f"  First {ref_name}: {bonded_sel_indices.size} bonds")
                            for sel_idx, dist_val in zip(bonded_sel_indices, bonded_distances):
                                logging.info(f"    -> {sel_name}[{sel_idx}]: {dist_val:.6f} nm")
                        else:
                            n_bonds = compute_bonds_vectorized(
                                ref_pos, sel_pos, box_frame, rcut_pair, calculate_lengths=False
                            )
                        frame_bonds_data.append(int(n_bonds))
                
                frame_bonds.append(frame_bonds_data)
                if args.BL:
                    frame_bond_lengths.append(frame_lengths_data)
                
                total_frames += 1
                
                # Skip frames based on skip parameter
                for _ in range(args.skip - 1):
                    read_gro_frame(f, skip=True)
        else:
            # Process specific frame range
            frame_idx = args.begin
            while frame_idx < args.end:
                coords, box_frame, _ = read_gro_frame(f, skip=False)
                if coords is None:
                    logging.warning(f"Reached EOF before frame {args.end}")
                    break
                
                # Process frame positions with MIC for molecules
                frame_data = process_frame_positions(args, groups, groups_mol, coords, box_frame)
                
                # Compute bonds for all pairs (explicit pairing: ref[i] with sel[i])
                frame_bonds_data = []
                frame_lengths_data = [] if args.BL else None
                
                for i, (ref_name, sel_name) in enumerate(zip(args.ref, args.sel)):
                    ref_pos = frame_data[('ref', i)]
                    sel_pos = frame_data[('sel', i)]
                    rcut_pair = args.rcut[i]
                    
                    # Debug: print bonds for first ref molecule/atom
                    debug_flag = args.debug and i == 0 and total_frames == 0
                    
                    if args.BL and args.dist:
                        if debug_flag:
                            result = compute_bonds_vectorized(
                                ref_pos, sel_pos, box_frame, rcut_pair, calculate_lengths=True, debug_first_ref=True, return_distances=True
                            )
                            n_bonds, avg_length, bonded_sel_indices, bonded_distances = result[:4]
                            logging.info(f"\n[DEBUG] Frame {total_frames}, Pair {ref_name}-{sel_name}:")
                            logging.info(f"  First {ref_name}: {bonded_sel_indices.size} bonds")
                            for sel_idx, dist_val in zip(bonded_sel_indices, bonded_distances):
                                logging.info(f"    -> {sel_name}[{sel_idx}]: {dist_val:.6f} nm")
                        else:
                            n_bonds, avg_length, pair_distances = compute_bonds_vectorized(
                                ref_pos, sel_pos, box_frame, rcut_pair, calculate_lengths=True, return_distances=True
                            )
                            bond_distances[i].extend(pair_distances)
                        frame_bonds_data.append(int(n_bonds))
                        frame_lengths_data.append(float(avg_length))
                    elif args.BL:
                        if debug_flag:
                            result = compute_bonds_vectorized(
                                ref_pos, sel_pos, box_frame, rcut_pair, calculate_lengths=True, debug_first_ref=True
                            )
                            n_bonds, avg_length, bonded_sel_indices, bonded_distances = result
                            logging.info(f"\n[DEBUG] Frame {total_frames}, Pair {ref_name}-{sel_name}:")
                            logging.info(f"  First {ref_name}: {bonded_sel_indices.size} bonds")
                            for sel_idx, dist_val in zip(bonded_sel_indices, bonded_distances):
                                logging.info(f"    -> {sel_name}[{sel_idx}]: {dist_val:.6f} nm")
                        else:
                            n_bonds, avg_length = compute_bonds_vectorized(
                                ref_pos, sel_pos, box_frame, rcut_pair, calculate_lengths=True
                            )
                        frame_bonds_data.append(int(n_bonds))
                        frame_lengths_data.append(float(avg_length))
                    elif args.dist:
                        if debug_flag:
                            result = compute_bonds_vectorized(
                                ref_pos, sel_pos, box_frame, rcut_pair, calculate_lengths=False, debug_first_ref=True, return_distances=True
                            )
                            n_bonds, bonded_sel_indices, bonded_distances = result[:3]
                            logging.info(f"\n[DEBUG] Frame {total_frames}, Pair {ref_name}-{sel_name}:")
                            logging.info(f"  First {ref_name}: {bonded_sel_indices.size} bonds")
                            for sel_idx, dist_val in zip(bonded_sel_indices, bonded_distances):
                                logging.info(f"    -> {sel_name}[{sel_idx}]: {dist_val:.6f} nm")
                        else:
                            n_bonds, pair_distances = compute_bonds_vectorized(
                                ref_pos, sel_pos, box_frame, rcut_pair, calculate_lengths=False, return_distances=True
                            )
                            bond_distances[i].extend(pair_distances)
                        frame_bonds_data.append(int(n_bonds))
                    else:
                        if debug_flag:
                            result = compute_bonds_vectorized(
                                ref_pos, sel_pos, box_frame, rcut_pair, calculate_lengths=False, debug_first_ref=True
                            )
                            n_bonds, bonded_sel_indices, bonded_distances = result
                            logging.info(f"\n[DEBUG] Frame {total_frames}, Pair {ref_name}-{sel_name}:")
                            logging.info(f"  First {ref_name}: {bonded_sel_indices.size} bonds")
                            for sel_idx, dist_val in zip(bonded_sel_indices, bonded_distances):
                                logging.info(f"    -> {sel_name}[{sel_idx}]: {dist_val:.6f} nm")
                        else:
                            n_bonds = compute_bonds_vectorized(
                                ref_pos, sel_pos, box_frame, rcut_pair, calculate_lengths=False
                            )
                        frame_bonds_data.append(int(n_bonds))
                
                frame_bonds.append(frame_bonds_data)
                if args.BL:
                    frame_bond_lengths.append(frame_lengths_data)
                
                total_frames += 1
                frame_idx += args.skip
                
                # Skip frames based on skip parameter
                for _ in range(args.skip - 1):
                    read_gro_frame(f, skip=True)
    
    # Convert to numpy arrays
    bond_counts = np.array(frame_bonds).T  # Shape: (n_pairs, n_frames)
    if args.BL:
        bond_lengths_avg = np.array(frame_bond_lengths).T  # Shape: (n_pairs, n_frames)
    else:
        bond_lengths_avg = None
    
    # Convert bond_distances lists to numpy arrays
    if args.dist:
        bond_distances = [np.array(d, dtype=np.float32) if d else None for d in bond_distances]
    
    return bond_counts, bond_lengths_avg, total_frames, bond_distances


def normalize_and_output(args, bond_counts, bond_lengths_avg, total_frames, ref_counts, sel_counts):
    """
    Normalize and write bond statistics to file(s).
    
    Outputs total number of bonds per frame for each pair.
    If --BL flag set, also outputs average bond length.
    """
    if total_frames == 0:
        logging.warning("No frames were processed. Exiting.")
        return
    
    # Create frame numbers for output
    frame_numbers = np.arange(total_frames)
    
    # Total bond counts (not normalized)
    output_data = [frame_numbers]
    header = "# frame"
    
    logging.info("--- Bond Count Output ---")
    logging.info(f"Total frames processed: {total_frames}")
    
    # Log rcut values
    logging.info("\n--- Cutoff Distance Configuration ---")
    if len(set(args.rcut)) == 1:
        logging.info(f"Cutoff distance (rcut): {args.rcut[0]:.3f} nm (all pairs)")
    else:
        logging.info(f"Cutoff distances (rcut) per pair:")
        for i, (ref_name, sel_name) in enumerate(zip(args.ref, args.sel)):
            logging.info(f"  {ref_name}-{sel_name}: {args.rcut[i]:.3f} nm")
    
    for i, (ref_name, sel_name) in enumerate(zip(args.ref, args.sel)):
        N_ref = int(ref_counts[i])
        N_sel = int(sel_counts[i])
        
        logging.info(f"Bond statistics for pair: {ref_name}-{sel_name}")
        logging.info(f"  N_ref ({ref_name}): {N_ref}")
        logging.info(f"  N_sel ({sel_name}): {N_sel}")
        
        # Total bonds per frame (not normalized)
        total_bonds = bond_counts[i, :]
        logging.info(f"  Mean total bonds per frame: {np.mean(total_bonds):.4f}")
        logging.info(f"  Min total bonds per frame: {np.min(total_bonds):.0f}")
        logging.info(f"  Max total bonds per frame: {np.max(total_bonds):.0f}")
        
        output_data.append(total_bonds)
        header += f" Total_bonds_{ref_name}-{sel_name}"
    
    # Write bond count output with header and footer comments
    with open(args.output_file, 'w') as f:
        # Write header
        f.write(header + "\n")
        # Write data
        np.savetxt(f, np.column_stack(output_data), fmt='%10.6f')
        # Write footer comments with counts
        f.write("# Total number of Ref:\n")
        f.write(f"# {', '.join([f'{name}: {int(count)}' for name, count in zip(args.ref, ref_counts)])}\n")
        f.write("# Total number of sel:\n")
        f.write(f"# {', '.join([f'{name}: {int(count)}' for name, count in zip(args.sel, sel_counts)])}\n")
    
    logging.info(f"Bond count data written to {args.output_file}")
    
    # Write bond length output if requested
    if args.BL and bond_lengths_avg is not None:
        output_data_bl = [frame_numbers]
        header_bl = "# frame"
        
        logging.info("\n--- Bond Length Statistics ---")
        
        for i, (ref_name, sel_name) in enumerate(zip(args.ref, args.sel)):
            logging.info(f"Bond length for pair: {ref_name}-{sel_name}")
            
            # Calculate average bond length (considering only frames with bonds)
            avg_lengths = bond_lengths_avg[i, :]
            
            # Filter out frames with no bonds for meaningful statistics
            valid_lengths = avg_lengths[avg_lengths > 0]
            
            if len(valid_lengths) > 0:
                logging.info(f"  Mean bond length: {np.mean(valid_lengths):.6f} nm")
                logging.info(f"  Min bond length: {np.min(valid_lengths):.6f} nm")
                logging.info(f"  Max bond length: {np.max(valid_lengths):.6f} nm")
                logging.info(f"  Frames with bonds: {len(valid_lengths)}/{total_frames}")
            else:
                logging.warning(f"  No bonds found for {ref_name}-{sel_name}")
            
            output_data_bl.append(avg_lengths)
            header_bl += f" BL_{ref_name}-{sel_name}(nm)"
        
        # Write bond length output
        bl_output_file = args.output_file.split('.')[0] + '_BL.' + args.output_file.split('.')[-1]
        np.savetxt(bl_output_file, np.column_stack(output_data_bl), header=header_bl, fmt='%10.6f')
        logging.info(f"Bond length data written to {bl_output_file}")


def generate_distance_distribution(args, bond_distances, ref_counts, sel_counts):
    """
    Generate and write bond distance distribution histograms.
    
    Args:
        args: command-line arguments (contains -dist N)
        bond_distances: (n_pairs, total_bonded_pairs) list of distance arrays
        ref_counts: reference particle/molecule counts
        sel_counts: selection particle/molecule counts
    """
    if args.dist is None or bond_distances is None:
        return
    
    logging.info("\n--- Bond Distance Distribution ---")
    
    n_pairs = len(args.ref)
    n_bins = args.dist
    
    # Create output file
    dist_output_file = args.output_file.split('.')[0] + '_dist.' + args.output_file.split('.')[-1]
    
    # Find the maximum rcut for x-axis
    max_rcut = max(args.rcut)
    
    # Create bin edges: from 0 to max(rcut)
    bin_edges = np.linspace(0, max_rcut, n_bins + 1)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
    
    output_data = [bin_centers]
    header = "# distance(nm)"
    
    for i, (ref_name, sel_name) in enumerate(zip(args.ref, args.sel)):
        rcut_pair = args.rcut[i]
        
        if bond_distances[i] is None or len(bond_distances[i]) == 0:
            hist = np.zeros(n_bins)
            logging.info(f"Distance distribution for pair {ref_name}-{sel_name}: No bonds found")
        else:
            # Create histogram with bin edges adjusted to the pair's rcut if needed
            hist, _ = np.histogram(bond_distances[i], bins=bin_edges)
            
            logging.info(f"Distance distribution for pair {ref_name}-{sel_name}:")
            logging.info(f"  Total bonded pairs: {len(bond_distances[i])}")
            logging.info(f"  Min distance: {np.min(bond_distances[i]):.6f} nm")
            logging.info(f"  Max distance: {np.max(bond_distances[i]):.6f} nm")
            logging.info(f"  Mean distance: {np.mean(bond_distances[i]):.6f} nm")
            logging.info(f"  Median distance: {np.median(bond_distances[i]):.6f} nm")
            logging.info(f"  Histogram bins: {n_bins}, range: 0.0-{max_rcut:.3f} nm")
        
        output_data.append(hist)
        header += f" hist_{ref_name}-{sel_name}"
    
    # Write distance distribution output
    with open(dist_output_file, 'w') as f:
        f.write(header + "\n")
        np.savetxt(f, np.column_stack(output_data), fmt='%12.6f')
        f.write("# Total number of Ref:\n")
        f.write(f"# {', '.join([f'{name}: {int(count)}' for name, count in zip(args.ref, ref_counts)])}\n")
        f.write("# Total number of sel:\n")
        f.write(f"# {', '.join([f'{name}: {int(count)}' for name, count in zip(args.sel, sel_counts)])}\n")
    
    logging.info(f"Distance distribution data written to {dist_output_file}")


def main():
    """Main function to run the bond analysis from GROMACS trajectory."""
    t_start = time.time()
    
    args = get_cli_args()
    need_mass_file = False

    if args.debug:
        pass  # Could enable more verbose logging

    # Handle rcut values and validate pair structure
    nref = len(args.ref)
    nsel = len(args.sel)
    
    # Validate that ref and sel have same length (explicit pairing by index)
    if nref != nsel:
        raise ValueError(
            f"--ref and --sel must have the same length for explicit pairing. "
            f"Got {nref} ref groups and {nsel} sel groups.")
    
    npairs = nref  # Number of pairs (ref[i] pairs with sel[i])
    
    # Validate and process rcut values
    if len(args.rcut) == 1:
        # Single rcut for all pairs
        args.rcut = args.rcut * npairs
    elif len(args.rcut) != npairs:
        raise ValueError(
            f"--rcut must have length 1 or {npairs} (number of pairs), "
            f"but got {len(args.rcut)}")
    
    nref = len(args.ref)
    if args.ref_mol is None:
        args.ref_mol = [False] * nref
    elif len(args.ref_mol) == 1:
        if args.ref_mol[0]:
            need_mass_file = True
        args.ref_mol = args.ref_mol * nref
    elif len(args.ref_mol) != nref:
        raise ValueError(
            f"--ref_mol must have length 1 or {nref}, "
            f"but got {len(args.ref_mol)}")
    else:
        need_mass_file = True

    nsel = len(args.sel)
    if args.sel_mol is None:
        args.sel_mol = [False] * nsel
    elif len(args.sel_mol) == 1:
        if args.sel_mol[0]:
            need_mass_file = True
        args.sel_mol = args.sel_mol * nsel
    elif len(args.sel_mol) != nsel:
        raise ValueError(
            f"--sel_mol must have length 1 or {nsel}, "
            f"but got {len(args.sel_mol)}")
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
    logging.info("BOND ANALYSIS FROM TRAJECTORY")
    logging.info("=" * 70)
    logging.info(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(f"Hostname: {hostname}")
    logging.info(f"Platform: {platform_info}")
    logging.info(f"Python: {python_version}")
    logging.info("")
    logging.info("Input Parameters:")
    for k, v in vars(args).items():
        logging.info(f"  {k}: {v}")
    
    logging.info("\n--- Initialization ---")
    logging.info("Parsing index file...")
    groups = parse_ndx(args.index_file)
    
    if args.debug:
        for group_name in args.ref:
            atom_indices = groups[group_name] + 1  # Convert back to 1-based
            logging.info(f"Reference group '{group_name}' has {len(groups[group_name])} atoms.")
            logging.info(f"  Atom indices: {atom_indices}")
        for group_name in args.sel:
            atom_indices = groups[group_name] + 1
            logging.info(f"Selection group '{group_name}' has {len(groups[group_name])} atoms.")
            logging.info(f"  Atom indices: {atom_indices}")
        
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
            if args.debug:
                for mol in groups_mol[ref_name][:3]:  # Show first 3 molecules
                    logging.info(f"  Mol {mol['mol_index']}: atoms {mol['atom_indices']}, masses {mol['atom_masses']}, total {mol['total_mass']:.3f}")
        else:
            ref_counts[i] = len(ref_indices)
            logging.info(f"No of atoms in ref group '{ref_name}': {ref_counts[i]:.0f}")
    
    for i, sel_name in enumerate(args.sel):
        sel_indices = groups[sel_name]
        if args.sel_mol[i]:
            sel_counts[i] = len(groups_mol[sel_name])
            logging.info(f"No of molecules in sel group '{sel_name}': {sel_counts[i]:.0f}")
            if args.debug:
                for mol in groups_mol[sel_name][:3]:
                    logging.info(f"  Mol {mol['mol_index']}: atoms {mol['atom_indices']}, masses {mol['atom_masses']}, total {mol['total_mass']:.3f}")
        else:
            sel_counts[i] = len(sel_indices)
            logging.info(f"No of atoms in sel group '{sel_name}': {sel_counts[i]:.0f}")
    
    logging.info("Reference group counts: " + ", ".join([f"{name}: {int(count)}" for name, count in zip(args.ref, ref_counts)]))
    logging.info("Selection group counts: " + ", ".join([f"{name}: {int(count)}" for name, count in zip(args.sel, sel_counts)]))
    
    if args.BL:
        logging.info("Bond length calculation enabled (--BL flag set)")

    logging.info("\n--- Calculation ---")
    logging.info("Computing bond statistics using vectorized distance calculations...")
    
    bond_counts, bond_lengths_avg, total_frames, bond_distances = calculate_bonds(args, groups, box, groups_mol)
    logging.info(f"Total frames processed: {total_frames}")
    
    logging.info("\n--- Output ---")
    normalize_and_output(args, bond_counts, bond_lengths_avg, total_frames, ref_counts, sel_counts)
    
    # Generate distance distribution if requested
    if args.dist:
        generate_distance_distribution(args, bond_distances, ref_counts, sel_counts)
    
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
    logging.info("Bond analysis complete!")
    logging.info("=" * 70)


if __name__ == '__main__':
    main()
