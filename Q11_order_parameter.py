#!/usr/bin/env python3
"""
Calculates the Q11 Steinhardt order parameter from a GROMACS trajectory.

The Q11 order parameter measures local crystalline ordering around each atom.
For each reference atom, neighbors within Rcut are identified and spherical 
harmonics are computed to characterize the local environment.

Reference: Steinhardt, Nelson, Ronchetti, PRL 47, 1297 (1981)
Q_6 = sqrt(4π/13 * Σ_{m=-6}^{6} |q_6^m|²)
"""
import argparse
import sys
import os
import logging
import time
import socket
import platform
import numpy as np
from scipy.special import sph_harm_y

def get_cli_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-s', '--structure_file', required=True, 
                        help='Structure file (.gro) to initialize the system.')
    parser.add_argument('-f', '--traj_file', required=True, 
                        help='Trajectory file (.gro format).')
    parser.add_argument('-n', '--index_file', required=True, 
                        help='Index file (.ndx).')
    parser.add_argument('-o', '--output_file', default='Q11.xvg', 
                        help='Output file name for Q11 results (default: Q11.xvg).')
    parser.add_argument('--ref', required=True, nargs='+', 
                        help='Name(s) of the reference group(s) in the index file.')
    parser.add_argument('--sel', required=True, nargs='+', 
                        help='Name(s) of the selection group(s) in the index file (neighbors).')
    parser.add_argument('--ref_mol', type=str2bool, nargs='*', default=None, 
                        help='Use center-of-mass (COM) for reference group(s).')
    parser.add_argument('--sel_mol', type=str2bool, nargs='*', default=None, 
                        help='Use center-of-mass (COM) for selection group(s).')
    parser.add_argument('--mass_file', type=str, default=None, 
                        help='File with atomic masses for COM calculation (required if using COM).')
    parser.add_argument('--rcut', type=float, nargs='*', default=None, 
                        help='Cutoff distance for nearest neighbors (nm). Can be: single value (for all pairs), multiple values (one per pair), or none (defaults to min(box)/2).')
    parser.add_argument('--begin', type=int, default=0, 
                        help='First frame to read from trajectory.')
    parser.add_argument('--end', type=int, default=-1, 
                        help='Last frame to read from trajectory. Use -1 for end.')
    parser.add_argument('--skip', type=int, default=1, 
                        help='Read every Nth frame.')
    parser.add_argument('--detail', action='store_true', 
                        help='Write detailed neighbor list output.')
    parser.add_argument('--debug', action='store_true', 
                        help='Enable debug mode.')
    
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
    """Parse a GROMACS index file (.ndx) and return a dictionary of atom groups."""
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
    """Parse a mass file with format: [group_name] atom_name = mass."""
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
            # Extract: molno(0:5), resname(5:10), atomname(10:15), x(20:28), y(28:36), z(36:44)
            atom_info.append((line[0:5].strip(), line[5:10].strip(), line[10:15].strip(), line[15:20].strip()))
            coords[i] = [float(line[20:28]), float(line[28:36]), float(line[36:44])]
        
        box_line = f.readline()
        box = np.array([float(x) for x in box_line.split()], dtype=np.float32)
        
        return coords, box, atom_info
    
    except (IOError, ValueError) as e:
        print(f"Error reading frame: {e}", file=sys.stderr)
        return None, None, None


def get_com(coords, groups_mol, mol_name):
    """
    Calculate center-of-mass for molecules in a group.
    
    Args:
        coords: (N_atoms, 3) coordinate array
        groups_mol: dict of molecular groups with mass info
        mol_name: name of the group to compute COM for
    
    Returns:
        (N_mols, 3) array of COM coordinates
    """
    n_mols = len(groups_mol[mol_name])
    com = np.zeros((n_mols, 3), dtype=np.float32)
    for mol_id, mol in enumerate(groups_mol[mol_name]):
        indices = np.array(mol['atom_indices'], dtype=int)
        coords_mol = coords[indices]
        masses = np.array(mol['atom_masses'], dtype=np.float32)
        com[mol_id] = np.sum(coords_mol * masses[:, np.newaxis], axis=0) / mol['total_mass']
    
    return com


def cartesian_to_spherical(delta):
    """
    Convert Cartesian coordinates to spherical (r, theta, phi).
    theta: angle from +z axis (0 to π)
    phi: azimuthal angle (0 to 2π)
    
    Args:
        delta: (3,) Cartesian vector
    
    Returns:
        (theta, phi) in radians
    """
    x, y, z = delta[0], delta[1], delta[2]
    r = np.sqrt(x**2 + y**2 + z**2)
    
    if r == 0:
        return 0, 0
    
    theta = np.arccos(z / r)  # angle from +z axis
    phi = np.arctan2(y, x)    # azimuthal angle
    
    return theta, phi


def compute_spherical_harmonics(theta, phi, l=11):
    """
    Compute spherical harmonics Y_l^m(theta, phi) for all m.
    
    Args:
        theta: zenith angle (0 to π)
        phi: azimuthal angle (0 to 2π)
        l: order (default 6)
    
    Returns:
        dict with keys -l to +l, values = Y_l^m
    """
    harmonics = {}
    for m in range(-l, l + 1):
        # scipy.special.sph_harm_y(m, l, phi, theta) - note order!
        harmonics[m] = sph_harm_y(m, l, phi, theta)
    
    return harmonics


def build_molecule_mapping(atom_info, atom_indices):
    """
    Build mapping from trajectory atom index to molecule number.
    Molecules are numbered sequentially (1, 2, 3, ...) based on residue.
    """
    atom_to_mol = {}  # {traj_atom_idx: molecule_number}
    mol_first_atom = {}  # {molecule_number: first_traj_atom_idx}
    
    # Get unique residues in this group in order
    # atom_info[i] = (residue_nr, resname, atomname, atom_nr)
    # Use atom_info[traj_idx][0] for residue number, NOT [3]
    unique_residues = []
    for traj_idx in atom_indices:
        res_nr = int(atom_info[traj_idx][0])  # FIXED: use [0], not [3]
        if res_nr not in unique_residues:
            unique_residues.append(res_nr)
    
    # Sort residues to ensure consistent numbering
    unique_residues_sorted = sorted(unique_residues)
    
    # Map each residue to a molecule number
    residue_to_mol = {res: mol_num for mol_num, res in enumerate(unique_residues_sorted, start=1)}
    
    # Map atoms to molecules and find first atom of each molecule
    for traj_idx in atom_indices:
        res_nr = int(atom_info[traj_idx][0])  # FIXED: use [0], not [3]
        mol_num = residue_to_mol[res_nr]
        atom_to_mol[traj_idx] = mol_num
        
        # Track first atom of molecule (lowest traj index)
        if mol_num not in mol_first_atom or traj_idx < mol_first_atom[mol_num]:
            mol_first_atom[mol_num] = traj_idx
    
    # Debug: Log molecule structure
    logging.debug(f"Molecule Mapping: {len(mol_first_atom)} unique molecules")
    for mol_id in sorted(list(mol_first_atom.keys())[:5]):  # First 5 molecules
        first_atom_traj_idx = mol_first_atom[mol_id]
        logging.debug(f"  Mol {mol_id}: first_atom_traj_idx={first_atom_traj_idx}, 1-indexed={first_atom_traj_idx+1}")
    
    return atom_to_mol, mol_first_atom


def compute_Q11_for_reference(neighbors_data, l=11):
    """
    Compute Q11 order parameter for a single reference atom.
    
    Args:
        neighbors_data: list of (delta_vector, distance) tuples for neighbors within rcut
        l: order (default 6)
    
    Returns:
        Q11 value (float), or np.nan if no neighbors
    """
    if len(neighbors_data) == 0:
        return np.nan
    
    # Initialize spherical harmonic accumulator
    q_lm = {}
    for m in range(-l, l + 1):
        q_lm[m] = 0.0 + 0.0j  # complex
    
    # Sum over all neighbors
    for delta, _ in neighbors_data:
        theta, phi = cartesian_to_spherical(delta)
        harmonics = compute_spherical_harmonics(theta, phi, l)
        for m in range(-l, l + 1):
            q_lm[m] += harmonics[m]
    
    # Compute Q_l
    # Q_l = sqrt(4π/(2l+1) * Σ_m |q_l^m|²)
    sum_sq = sum(np.abs(q_lm[m])**2 for m in range(-l, l + 1))
    Q_l = np.sqrt(4.0 * np.pi / (2 * l + 1) * sum_sq)
    
    return Q_l

def calculate_q6_trajectory(args, groups, is_self_pair, groups_mol, box_init):
    """
    Calculate Q11 order parameter for each frame (per-pair tracking).
    """
    npairs = len(args.ref)
    q6_timeseries_by_pair = {i: [] for i in range(npairs)}
    nn_timeseries_by_pair = {i: [] for i in range(npairs)}
    q6_per_ref = {}
    neighbor_stats = {'n_neighbors': [], 'n_ref_atoms': []}
    all_neighbor_details = []
    
    frame_number = 0
    with open(args.traj_file, 'r') as f:
        for _ in range(args.begin):
            read_gro_frame(f, skip=True)
        
        if args.end == -1:
            while True:
                coords, box, atom_info = read_gro_frame(f, skip=False)
                if coords is None:
                    break
                
                if frame_number % max(1, args.skip) == 0:
                    atom_info = np.array(atom_info)
                    q6_by_pair, nn_by_pair, neighbor_details = process_frame_q6(args, groups, is_self_pair, groups_mol, 
                                                coords, box, atom_info, q6_per_ref, frame_number)
                    for pair_idx in range(npairs):
                        q6_timeseries_by_pair[pair_idx].append(q6_by_pair[pair_idx])
                        nn_timeseries_by_pair[pair_idx].append(nn_by_pair[pair_idx])
                    all_neighbor_details.extend(neighbor_details)
                    neighbor_stats['n_ref_atoms'].append(len(q6_per_ref))
                
                frame_number += 1
                for _ in range(args.skip - 1):
                    read_gro_frame(f, skip=True)
                    frame_number += 1
        else:
            for frame_idx in range(args.begin, args.end, args.skip):
                coords, box, atom_info = read_gro_frame(f, skip=False)
                if coords is None:
                    logging.warning(f"Reached EOF before frame {args.end}")
                    break
                
                atom_info = np.array(atom_info)
                q6_by_pair, nn_by_pair, neighbor_details = process_frame_q6(args, groups, is_self_pair, groups_mol, 
                                            coords, box, atom_info, q6_per_ref, frame_number)
                for pair_idx in range(npairs):
                    q6_timeseries_by_pair[pair_idx].append(q6_by_pair[pair_idx])
                    nn_timeseries_by_pair[pair_idx].append(nn_by_pair[pair_idx])
                all_neighbor_details.extend(neighbor_details)
                neighbor_stats['n_ref_atoms'].append(len(q6_per_ref))
                frame_number += 1
    
    for pair_idx in range(npairs):
        q6_timeseries_by_pair[pair_idx] = np.array(q6_timeseries_by_pair[pair_idx])
        nn_timeseries_by_pair[pair_idx] = np.array(nn_timeseries_by_pair[pair_idx])
    
    return q6_timeseries_by_pair, nn_timeseries_by_pair, q6_per_ref, neighbor_stats, all_neighbor_details



def process_frame_q6(args, groups, is_self_pair, groups_mol, coords, box, atom_info, 
                     q6_per_ref, frame_idx):
    """
    Process a single frame to compute Q11 for all reference atoms.
    Iterate over molecules (not atoms).
    Report each molecule (1, 2, 3, ...) with its neighbors.
    """
    npairs = len(args.ref)
    q6_by_pair = {}
    nn_by_pair = {}
    neighbor_details = []
    
    for pair_idx, (ref_name, sel_name) in enumerate(zip(args.ref, args.sel)):
        ref_indices = groups[ref_name]
        sel_indices = groups[sel_name]
        
        # Build per-pair molecule mappings
        ref_atom_to_mol, ref_mol_first_atom = build_molecule_mapping(atom_info, ref_indices)
        sel_atom_to_mol, sel_mol_first_atom = build_molecule_mapping(atom_info, sel_indices)
        
        if args.ref_mol[pair_idx]:
            ref_pos = get_com(coords, groups_mol, ref_name)
        else:
            ref_pos = coords[ref_indices]
        
        if args.sel_mol[pair_idx]:
            sel_pos = get_com(coords, groups_mol, sel_name)
        else:
            sel_pos = coords[sel_indices]
        
        pair_q6_values = []
        pair_nn_values = []
        
        # Iterate over ALL reference molecules (not atoms)
        sorted_all_ref_mols = sorted(ref_mol_first_atom.keys())  # All molecule IDs
        sorted_all_sel_mols = sorted(sel_mol_first_atom.keys())  # All molecule IDs
        
        # Create mapping from molecule ID to COM position index
        # (only needed when using COM positions)
        ref_mol_to_comidx = {mol_id: idx for idx, mol_id in enumerate(sorted_all_ref_mols)}
        sel_mol_to_comidx = {mol_id: idx for idx, mol_id in enumerate(sorted_all_sel_mols)}
        
        # For atom-based coordinates, create atom-to-posidx mapping
        if not args.ref_mol[pair_idx]:
            ref_atom_to_posidx = {int(ref_indices[i]): i for i in range(len(ref_indices))}
        if not args.sel_mol[pair_idx]:
            sel_atom_to_posidx = {int(sel_indices[i]): i for i in range(len(sel_indices))}
        
        for ref_mol_id_global in sorted_all_ref_mols:
            ref_mol_id_in_pair = sorted_all_ref_mols.index(ref_mol_id_global) + 1  # 1, 2, 3, ...
            
            # Get reference position (COM or atom)
            if args.ref_mol[pair_idx]:
                # Using COM: use molecule index directly
                ref_com_idx = ref_mol_to_comidx[ref_mol_id_global]
                r_pos = ref_pos[ref_com_idx]
                ref_atoms_to_process = [None]  # Single position per molecule
            else:
                # Using atoms: find all atoms in this molecule
                ref_atoms_to_process = [atom_idx for atom_idx, mol_id in ref_atom_to_mol.items() 
                                       if mol_id == ref_mol_id_global]
            
            neighbors = []
            neighbor_molecules = {}  # sel_mol_id_global -> (info)
            
            # For each atom (or pseudo-atom for COM) in reference molecule
            for ref_atom_idx in ref_atoms_to_process:
                if args.ref_mol[pair_idx]:
                    # Already have COM position
                    r_pos_list = [r_pos]
                else:
                    # Get atom position
                    if ref_atom_idx not in ref_atom_to_posidx:
                        continue
                    ref_pos_idx = ref_atom_to_posidx[ref_atom_idx]
                    r_pos_list = [ref_pos[ref_pos_idx]]
                
                # Check all selection molecules
                for sel_mol_id_global in sorted_all_sel_mols:
                    # Skip self-pair (only for molecules, not atoms within)
                    if is_self_pair[pair_idx] and ref_mol_id_global == sel_mol_id_global:
                        continue
                    
                    # Get selection position (COM or atom)
                    if args.sel_mol[pair_idx]:
                        # Using COM: use molecule index directly
                        sel_com_idx = sel_mol_to_comidx[sel_mol_id_global]
                        s_pos_list = [sel_pos[sel_com_idx]]
                    else:
                        # Using atoms: find all atoms in this molecule
                        sel_atoms_in_mol = [atom_idx for atom_idx, mol_id in sel_atom_to_mol.items() 
                                           if mol_id == sel_mol_id_global]
                        s_pos_list = []
                        for sel_atom_idx in sel_atoms_in_mol:
                            if sel_atom_idx not in sel_atom_to_posidx:
                                continue
                            sel_pos_idx = sel_atom_to_posidx[sel_atom_idx]
                            s_pos_list.append(sel_pos[sel_pos_idx])
                    
                    # Compute distances to all selection positions
                    for s_pos in s_pos_list:
                        for r_pos_actual in r_pos_list:
                            delta = s_pos - r_pos_actual
                            delta -= box[:3] * np.round(delta / box[:3])
                            dist = np.linalg.norm(delta)
                            
                            if dist < args.rcut[pair_idx]:
                                neighbors.append((delta, dist))
                                
                                sel_mol_id_in_pair = sorted_all_sel_mols.index(sel_mol_id_global) + 1
                                sel_atomID = sel_mol_first_atom[sel_mol_id_global] + 1  # 1-indexed trajectory atom ID
                                
                                # Track closest distance for each molecule pair
                                if sel_mol_id_global not in neighbor_molecules:
                                    neighbor_molecules[sel_mol_id_global] = {
                                        'dist': dist, 
                                        'delta': delta,
                                        'sel_mol_id_in_pair': sel_mol_id_in_pair,
                                        'sel_atomID': sel_atomID
                                    }
                                elif dist < neighbor_molecules[sel_mol_id_global]['dist']:
                                    neighbor_molecules[sel_mol_id_global] = {
                                        'dist': dist, 
                                        'delta': delta,
                                        'sel_mol_id_in_pair': sel_mol_id_in_pair,
                                        'sel_atomID': sel_atomID
                                    }
            
            # Compute Q11 for this reference molecule
            q6 = compute_Q11_for_reference(neighbors, l=11)
            
            global_ref_id = len(ref_indices) * pair_idx + ref_mol_id_in_pair - 1
            if global_ref_id not in q6_per_ref:
                q6_per_ref[global_ref_id] = []
            q6_per_ref[global_ref_id].append((q6, len(neighbors)))
            
            pair_q6_values.append(q6)
            pair_nn_values.append(len(neighbors))
            
            # Report this molecule with its neighbors
            if neighbor_molecules:  # Only report if has neighbors
                ref_atomID = ref_mol_first_atom[ref_mol_id_global] + 1  # 1-indexed trajectory atom ID
                neighbors_list = []
                for sel_mol_id_global in sorted(neighbor_molecules.keys()):
                    info = neighbor_molecules[sel_mol_id_global]
                    neighbors_list.append((info['sel_mol_id_in_pair'], info['sel_atomID'], info['dist'], info['delta']))
                
                neighbor_details.append({
                    'frame': frame_idx,
                    'pair_idx': pair_idx,
                    'ref_name': ref_name,
                    'ref_id': ref_mol_id_in_pair,
                    'ref_atom_id': ref_atomID,
                    'sel_name': sel_name,
                    'neighbors': neighbors_list
                })
        
        pair_q6_array = np.array(pair_q6_values)
        q6_by_pair[pair_idx] = np.nanmean(pair_q6_array)
        nn_by_pair[pair_idx] = np.mean(pair_nn_values) if pair_nn_values else 0.0
    
    return q6_by_pair, nn_by_pair, neighbor_details


def write_q6_output(args, q6_timeseries_by_pair, nn_timeseries_by_pair, all_neighbor_details):
    """Write Q11 results to output file - multi-column format + enhanced detail."""
    npairs = len(q6_timeseries_by_pair)
    
    if npairs == 0 or len(q6_timeseries_by_pair[0]) == 0:
        logging.warning("No frames processed. No output written.")
        return
    
    base_name = os.path.splitext(args.output_file)[0]
    
    # Write Q11 file with multiple columns (one per pair)
    q6_header = "# Time(ps)" + "".join([f"    {args.ref[i]}-{args.sel[i]}" for i in range(npairs)])
    q6_output_data = []
    for frame_idx in range(len(q6_timeseries_by_pair[0])):
        time_ps = frame_idx * args.skip
        row = [time_ps]
        for pair_idx in range(npairs):
            row.append(q6_timeseries_by_pair[pair_idx][frame_idx])
        q6_output_data.append(row)
    
    q6_array = np.array(q6_output_data)
    fmt_str = f"{'%10.3f'} {' %10.6f' * npairs}"
    np.savetxt(args.output_file, q6_array, header=q6_header, fmt=fmt_str, comments='')
    logging.info(f"Q11 time series written to {args.output_file}")
    
    logging.info("\n--- Q11 Statistics ---")
    for pair_idx in range(npairs):
        pair_name = f"{args.ref[pair_idx]}-{args.sel[pair_idx]}"
        q6_data = q6_timeseries_by_pair[pair_idx]
        nn_data = nn_timeseries_by_pair[pair_idx]
        
        logging.info(f"\n{pair_name}:")
        logging.info(f"  Mean Q11:     {np.nanmean(q6_data):.6f}")
        logging.info(f"  Std Dev Q11:  {np.nanstd(q6_data):.6f}")
        logging.info(f"  Min Q11:      {np.nanmin(q6_data):.6f}")
        logging.info(f"  Max Q11:      {np.nanmax(q6_data):.6f}")
        logging.info(f"  Mean Avg_NN: {np.nanmean(nn_data):.2f}")
    
    # Write NN file
    nn_file = f"{base_name}_NN.xvg"
    nn_header = "# Time(ps)" + "".join([f"    {args.ref[i]}-{args.sel[i]}" for i in range(npairs)])
    nn_output_data = []
    for frame_idx in range(len(q6_timeseries_by_pair[0])):
        time_ps = frame_idx * args.skip
        row = [time_ps]
        for pair_idx in range(npairs):
            row.append(nn_timeseries_by_pair[pair_idx][frame_idx])
        nn_output_data.append(row)
    
    nn_array = np.array(nn_output_data)
    fmt_str = f"{'%10.3f'} {' %10.2f' * npairs}"
    np.savetxt(nn_file, nn_array, header=nn_header, fmt=fmt_str, comments='')
    logging.info(f"\nAverage NN time series written to {nn_file}")
    
    if args.detail:
        detail_file = os.path.splitext(args.output_file)[0] + '_detail.txt'
        with open(detail_file, 'w') as f:
            f.write("# Frame  Time(ps)  Pair_id  Ref_id  Ref_AtomID  Ref_Name  NN_id  Sel_id  Sel_AtomID  Sel_Name  Distance(nm)  Vec_X  Vec_Y  Vec_Z  Total_NN\n")
            for detail in all_neighbor_details:
                frame = detail['frame']
                time_ps = frame * args.skip
                pair_id = detail['pair_idx'] + 1
                ref_name = detail['ref_name']
                ref_id = detail['ref_id']
                ref_atom_id = detail['ref_atom_id']
                sel_name = detail['sel_name']
                total_nn = len(detail['neighbors'])
                for nn_idx, (sel_id, sel_atom_id, distance, vector) in enumerate(detail['neighbors'], start=1):
                    f.write(f"{frame:6d} {time_ps:10.3f} {pair_id:8d} {ref_id:8d} {ref_atom_id:11d} {ref_name:8s} {nn_idx:6d} {sel_id:8d} {sel_atom_id:11d} {sel_name:8s} {distance:13.6f} {vector[0]:10.6f} {vector[1]:10.6f} {vector[2]:10.6f} {total_nn:9d}\n")
        logging.info(f"\nDetailed information written to {detail_file}")
def main():
    """Main function to run Q11 order parameter calculation."""
    # Record start time
    start_time = time.time()
    
    args = get_cli_args()
    
    # Validate and expand ref_mol/sel_mol arguments
    nref = len(args.ref)
    if args.ref_mol is None:
        args.ref_mol = [False] * nref
    elif len(args.ref_mol) == 1:
        args.ref_mol = args.ref_mol * nref
    elif len(args.ref_mol) != nref:
        raise ValueError(f"--ref_mol must have length 1 or {nref}, got {len(args.ref_mol)}")
    
    nsel = len(args.sel)
    if nsel != nref:
        raise ValueError(f"Number of selection groups ({nsel}) does not match number of reference groups ({nref}).")

    if args.sel_mol is None:
        args.sel_mol = [False] * nsel
    elif len(args.sel_mol) == 1:
        args.sel_mol = args.sel_mol * nsel
    elif len(args.sel_mol) != nsel:
        raise ValueError(f"--sel_mol must have length 1 or {nsel}, got {len(args.sel_mol)}")
    
    npair=nref
    # Check if mass file is needed
    need_mass = any(args.ref_mol) or any(args.sel_mol)
    if need_mass and args.mass_file is None:
        raise ValueError("Mass file required for COM calculations but not provided.")
    
    # Setup logging
    log_file = os.path.splitext(args.output_file)[0] + '.log'
    logging.basicConfig(
        level=logging.INFO,
        format='%(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    logging.info("=" * 70)
    logging.info("Q11 STEINHARDT ORDER PARAMETER CALCULATION")
    logging.info("=" * 70)
    logging.info("Command-line parameters:")
    for k, v in vars(args).items():
        if k not in ['ref', 'sel', 'ref_mol', 'sel_mol', 'rcut']:
            logging.info(f"  {k}: {v}")
    
    # Parse index file
    logging.info("\n--- Initialization ---")
    logging.info("Parsing index file...")
    groups = parse_ndx(args.index_file)
    
    # Log atom counts for each group
    logging.info("\nGroup Summary:")
    for i, ref_name in enumerate(args.ref):
        n_ref_atoms = len(groups[ref_name])
        sel_name = args.sel[i]
        n_sel_atoms = len(groups[sel_name])
        logging.info(f"  Pair {i+1}: {ref_name} ({n_ref_atoms} atoms) <-> {sel_name} ({n_sel_atoms} atoms)")
    
    if args.debug:
        for group_name in args.ref + args.sel:
            if group_name in groups:
                logging.info(f"  DEBUG: Group '{group_name}' atom IDs (1-indexed): {groups[group_name] + 1}")
    
    # Parse mass file if needed
    groups_mol = None
    if need_mass:
        logging.info("\nParsing mass file...")
        masses = parse_mass_file(args.mass_file)
        logging.info(f"Found {len(masses)} molecular groups in mass file")
    
    # Read structure file to initialize molecular groups
    logging.info("Reading structure file...")
    logging.info(f"  File: {args.structure_file}")
    with open(args.structure_file, 'r') as f:
        coords_init, box_init, atom_info_init = read_gro_frame(f)
    
    if coords_init is None:
        logging.error("Failed to read structure file.")
        sys.exit(1)
    
    logging.info(f"  Total atoms: {len(coords_init)}")
    logging.info(f"  Box dimensions: {box_init[0]:.4f} x {box_init[1]:.4f} x {box_init[2]:.4f} nm")
    
    # Build molecular groups for COM calculations
    if need_mass:
        atom_info_init = np.array(atom_info_init)
        groups_mol = {}
        for group_name in args.ref + args.sel:
            if group_name not in groups:
                continue
            atom_indices = groups[group_name]
            
            mols_in_group = {}
            for idx in atom_indices:
                mol_no = atom_info_init[idx][0]
                atom_name = atom_info_init[idx][2]
                
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
                    logging.warning(f"Mass not found for {group_name}/{atom_name}, using 0.0")
                    mass = 0.0
                
                mols_in_group[mol_no]["atom_masses"].append(mass)
                mols_in_group[mol_no]["total_mass"] = np.sum(mols_in_group[mol_no]["atom_masses"])
            
            groups_mol[group_name] = list(mols_in_group.values())
            n_mols = len(groups_mol[group_name])
            logging.info(f"\nGroup '{group_name}': {n_mols} molecules")
            
            # Log molecule details
            if n_mols > 0:
                mol_data = groups_mol[group_name][0]
                n_atoms_per_mol = len(mol_data['atom_indices'])
                total_mass = mol_data['total_mass']
                logging.info(f"  Atoms per molecule: {n_atoms_per_mol}")
                logging.info(f"  Atom types: {', '.join(set(mol_data['atom_names']))}")
                logging.info(f"  Mass per molecule: {total_mass:.4f} u")
                if args.debug:
                    logging.info(f"  DEBUG: First molecule atoms: {mol_data['atom_names']} with masses {mol_data['atom_masses']}")
    
    # Set rcut if not provided or expand if needed
    if args.rcut is None or len(args.rcut) == 0:
        # Use min(box)/2 for all pairs
        args.rcut = [np.min(box_init[:3]) / 2.0] * npair
        logging.info(f"Rcut not provided, using min(box)/2 = {args.rcut[0]:.4f} nm for all {npair} pair(s)")
    elif len(args.rcut) == 1:
        # Use single value for all pairs
        args.rcut = args.rcut * npair
        logging.info(f"Using single rcut = {args.rcut[0]:.4f} nm for all {npair} pair(s)")
    elif len(args.rcut) == npair:
        # One value per pair
        logging.info(f"Using pair-specific rcut values:")
        for i, (ref_name, sel_name) in enumerate(zip(args.ref, args.sel)):
            logging.info(f"  Pair {i+1} ({ref_name}-{sel_name}): rcut = {args.rcut[i]:.4f} nm")
    else:
        raise ValueError(f"Number of rcut values ({len(args.rcut)}) does not match number of pairs ({npair}). Expected 1, {npair}, or none.")
    
    # Convert to list if it's still a single value after the loop
    if not isinstance(args.rcut, list):
        args.rcut = list(args.rcut)
    
    # Determine self-pair status (1-to-1 pairing)
    is_self_pair = np.zeros(npair, dtype=bool)
    logging.info("\nPair Configuration:")
    for i, (ref_name, sel_name) in enumerate(zip(args.ref, args.sel)):
        if ref_name == sel_name:
            is_self_pair[i] = True
        self_pair_str = "YES (self-comparison)" if is_self_pair[i] else "NO (cross-pair)"
        logging.info(f"  Pair {i+1}: {ref_name} <-> {sel_name}")
        logging.info(f"    Self-pair: {self_pair_str}")
        logging.info(f"    Cutoff: {args.rcut[i]:.4f} nm")
        logging.info(f"    COM mode: ref={args.ref_mol[i]}, sel={args.sel_mol[i]}")
    
    # Calculate Q11 for trajectory
    logging.info("\n--- Calculation ---")
    logging.info("Processing trajectory...")
    q6_timeseries_by_pair, nn_timeseries_by_pair, q6_per_ref, neighbor_stats, all_neighbor_details = calculate_q6_trajectory(
        args, groups, is_self_pair, groups_mol, box_init
    )
    
    # Write output
    logging.info("\n--- Output ---")
    write_q6_output(args, q6_timeseries_by_pair, nn_timeseries_by_pair, all_neighbor_details)
    
    logging.info("\n" + "=" * 70)
    logging.info("CALCULATION COMPLETED")
    logging.info("=" * 70)
    
    # Print computation details
    end_time = time.time()
    elapsed_time = end_time - start_time
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = elapsed_time % 60
    
    logging.info("\n--- Computation Details ---")
    logging.info(f"Start time:    {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start_time))}")
    logging.info(f"End time:      {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end_time))}")
    logging.info(f"Elapsed time:  {hours}h {minutes}m {seconds:.2f}s")
    logging.info(f"Machine name:  {socket.gethostname()}")
    logging.info(f"OS:            {platform.system()} {platform.release()}")
    logging.info(f"Python:        {platform.python_version()}")
    logging.info(f"Processor:     {platform.processor()}")
    
    


if __name__ == '__main__':
    main()
