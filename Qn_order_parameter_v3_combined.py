#!/usr/bin/env python3
"""
Qn STEINHARDT ORDER PARAMETER - V3 COMBINED
Calculates the Qn Steinhardt order parameter from a GROMACS trajectory.

OPTIMIZATION: Vectorized operations + Frame-level parallelization
- Combines cdist (vectorized distances) with multiprocessing.Pool
- Each worker uses vectorized spherical harmonics computation
- Maximum performance: ~2-3x from vectorization, ~3-4x from parallelization = 6-12x total
- Excellent for multiple pairs - each worker handles all pairs for a frame simultaneously
- Order parameter n can be specified per pair via --order.

Reference: Steinhardt, Nelson, Ronchetti, PRL 47, 1297 (1981)
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
from scipy.spatial.distance import cdist
from multiprocessing import Pool
import os

def get_cli_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-s', '--structure_file', required=True, help='Structure file (.gro)')
    parser.add_argument('-f', '--traj_file', required=True, help='Trajectory file (.gro)')
    parser.add_argument('-n', '--index_file', required=True, help='Index file (.ndx)')
    parser.add_argument('-o', '--output_file', default='Qn.xvg', help='Output file')
    parser.add_argument('--order', type=int, nargs='+', default=None, help='Order parameter l for each pair')
    parser.add_argument('--ref', required=True, nargs='+', help='Reference group(s)')
    parser.add_argument('--sel', required=True, nargs='+', help='Selection group(s)')
    parser.add_argument('--ref_mol', type=str2bool, nargs='*', default=None, help='Use COM for ref')
    parser.add_argument('--sel_mol', type=str2bool, nargs='*', default=None, help='Use COM for sel')
    parser.add_argument('--mass_file', type=str, default=None, help='Mass file for COM')
    parser.add_argument('--rcut', type=float, nargs='*', default=None, help='Cutoff distance (nm)')
    parser.add_argument('--begin', type=int, default=0, help='First frame')
    parser.add_argument('--end', type=int, default=-1, help='Last frame')
    parser.add_argument('--skip', type=int, default=1, help='Read every Nth frame')
    parser.add_argument('--detail', action='store_true', help='Write detailed output')
    parser.add_argument('--debug', action='store_true', help='Debug mode')
    parser.add_argument('-j', '--nproc', type=int, default=None, help='Number of processes')
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
        raise argparse.ArgumentTypeError('Boolean value expected')

def parse_ndx(ndx_file):
    """Parse GROMACS index file."""
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
    for name, indices in groups.items():
        groups[name] = np.array(indices, dtype=int) - 1
    return groups

def parse_mass_file(mass_file):
    """Parse mass file."""
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
    """Read frame from GRO file."""
    title = f.readline()
    if not title:
        return None, None, None
    try:
        num_atoms_str = f.readline()
        if not num_atoms_str:
            return None, None, None
        num_atoms = int(num_atoms_str.strip())
        if skip:
            for _ in range(num_atoms + 1):
                f.readline()
            return None, None, None
        coords = np.zeros((num_atoms, 3), dtype=np.float32)
        atom_info = []
        for i in range(num_atoms):
            line = f.readline()
            atom_info.append((line[0:5].strip(), line[5:10].strip(), line[10:15].strip(), line[15:20].strip()))
            coords[i] = [float(line[20:28]), float(line[28:36]), float(line[36:44])]
        box_line = f.readline()
        box = np.array([float(x) for x in box_line.split()], dtype=np.float32)
        return coords, box, atom_info
    except (IOError, ValueError) as e:
        print(f"Error reading frame: {e}", file=sys.stderr)
        return None, None, None

def get_com(coords, groups_mol, mol_name):
    """Calculate center-of-mass."""
    n_mols = len(groups_mol[mol_name])
    com = np.zeros((n_mols, 3), dtype=np.float32)
    for mol_id, mol in enumerate(groups_mol[mol_name]):
        indices = np.array(mol['atom_indices'], dtype=int)
        coords_mol = coords[indices]
        masses = np.array(mol['atom_masses'], dtype=np.float32)
        com[mol_id] = np.sum(coords_mol * masses[:, np.newaxis], axis=0) / mol['total_mass']
    return com

def cartesian_to_spherical_vectorized(deltas):
    """Vectorized Cartesian to spherical."""
    x = deltas[:, 0]
    y = deltas[:, 1]
    z = deltas[:, 2]
    r = np.sqrt(x**2 + y**2 + z**2)
    r_safe = np.where(r == 0, 1.0, r)
    theta = np.arccos(z / r_safe)
    phi = np.arctan2(y, x)
    return theta, phi

def compute_q6_vectorized(deltas, l=6):
    """Compute Q6 using vectorized harmonics."""
    if len(deltas) == 0:
        return np.nan
    theta, phi = cartesian_to_spherical_vectorized(deltas)
    q_lm = np.zeros((2*l + 1,), dtype=np.complex128)
    for m_idx, m in enumerate(range(-l, l + 1)):
        harmonics = sph_harm_y(m, l, phi, theta)
        q_lm[m_idx] = np.sum(harmonics)
    sum_sq = np.sum(np.abs(q_lm)**2)
    Q_l = np.sqrt(4.0 * np.pi / (2 * l + 1) * sum_sq)
    return Q_l

def apply_mic_vectorized(deltas, box):
    """Apply MIC to batch."""
    return deltas - box[:3] * np.round(deltas / box[:3])

def process_frame_q6_combined(frame_data_tuple):
    """
    Worker function combining vectorized operations + parallelization.
    Reports results by molecule.
    
    Args:
        frame_data_tuple: (coords, box, atom_info, args_dict, groups_dict, is_self_pair, groups_mol, frame_idx)
    """
    coords, box, atom_info, args_dict, groups_dict, is_self_pair, groups_mol, frame_idx = frame_data_tuple
    
    neighbor_details = []
    q6_per_ref_frame = {}
    q6_by_pair = {}  # Track Q6 separately for each pair
    nn_by_pair = {}  # Track NN separately for each pair
    
    atom_info = np.array(atom_info)
    
    for i, (ref_name, sel_name) in enumerate(zip(args_dict['ref'], args_dict['sel'])):
        ref_indices = groups_dict[ref_name]
        sel_indices = groups_dict[sel_name]
        
        # Build molecule mappings
        ref_atom_to_mol = {}
        ref_mol_first_atom = {}
        unique_ref_residues = []
        for traj_idx in ref_indices:
            res_nr = int(atom_info[traj_idx][0])
            if res_nr not in unique_ref_residues:
                unique_ref_residues.append(res_nr)
        # Sort residues to match V0 behavior
        unique_ref_residues_sorted = sorted(unique_ref_residues)
        residue_to_mol_ref = {res: mol_num for mol_num, res in enumerate(unique_ref_residues_sorted, start=1)}
        for traj_idx in ref_indices:
            res_nr = int(atom_info[traj_idx][0])
            mol_num = residue_to_mol_ref[res_nr]
            ref_atom_to_mol[traj_idx] = mol_num
            if mol_num not in ref_mol_first_atom or traj_idx < ref_mol_first_atom[mol_num]:
                ref_mol_first_atom[mol_num] = traj_idx
        
        # Same for selection
        sel_atom_to_mol = {}
        sel_mol_first_atom = {}
        unique_sel_residues = []
        for traj_idx in sel_indices:
            res_nr = int(atom_info[traj_idx][0])
            if res_nr not in unique_sel_residues:
                unique_sel_residues.append(res_nr)
        # Sort residues to match V0 behavior
        unique_sel_residues_sorted = sorted(unique_sel_residues)
        residue_to_mol_sel = {res: mol_num for mol_num, res in enumerate(unique_sel_residues_sorted, start=1)}
        for traj_idx in sel_indices:
            res_nr = int(atom_info[traj_idx][0])
            mol_num = residue_to_mol_sel[res_nr]
            sel_atom_to_mol[traj_idx] = mol_num
            if mol_num not in sel_mol_first_atom or traj_idx < sel_mol_first_atom[mol_num]:
                sel_mol_first_atom[mol_num] = traj_idx
        
        if args_dict['ref_mol'][i]:
            ref_pos = get_com(coords, groups_mol, ref_name)
        else:
            ref_pos = coords[ref_indices]
        
        if args_dict['sel_mol'][i]:
            sel_pos = get_com(coords, groups_mol, sel_name)
        else:
            sel_pos = coords[sel_indices]
        
        # VECTORIZED: Compute all distances at once
        deltas_full = sel_pos[np.newaxis, :, :] - ref_pos[:, np.newaxis, :]  # (n_ref, n_sel, 3)
        deltas_full = apply_mic_vectorized(deltas_full.reshape(-1, 3), box).reshape(deltas_full.shape)
        distance_matrix = np.linalg.norm(deltas_full, axis=2)
        
        pair_q6_values = []  # Collect Q6 for this pair only
        pair_nn_values = []  # Collect NN for this pair only
        mol_neighbors = {}  # Track neighbors for each ref molecule: ref_mol_id -> {sel_mol_id -> info}
        
        # Get sorted molecule lists upfront (same as V0)
        sorted_all_ref_mols = sorted(ref_mol_first_atom.keys())
        sorted_all_sel_mols = sorted(sel_mol_first_atom.keys())
        
        # Create mapping: position index -> mol_id (for COM case)
        ref_posidx_to_molid = {idx: mol_id for idx, mol_id in enumerate(sorted_all_ref_mols)}
        sel_posidx_to_molid = {idx: mol_id for idx, mol_id in enumerate(sorted_all_sel_mols)}
        
        # Process each reference atom/molecule
        for ref_atom_id in range(len(ref_pos)):
            neighbor_mask = distance_matrix[ref_atom_id] < args_dict['rcut'][i]
            neighbor_indices = np.where(neighbor_mask)[0]
            
            if is_self_pair[i]:
                neighbor_indices = neighbor_indices[neighbor_indices != ref_atom_id]
            
            # Get reference molecule information
            if args_dict['ref_mol'][i]:
                # Using COM: ref_atom_id is position index, map to mol_id
                ref_mol_id = ref_posidx_to_molid[ref_atom_id]
            else:
                # Using atoms: get molecule from atom index
                ref_atom_traj_id = ref_indices[ref_atom_id]
                ref_mol_id = ref_atom_to_mol[ref_atom_traj_id]
            
            # Filter neighbors: exclude atoms/molecules from same molecule if self-pair
            if is_self_pair[i]:
                filtered_indices = []
                for idx in neighbor_indices:
                    if args_dict['sel_mol'][i]:
                        # Using COM: idx is position index, map to mol_id
                        sel_mol_id = sel_posidx_to_molid[idx]
                    else:
                        # Using atoms: get molecule from atom index
                        sel_atom_traj_id = sel_indices[idx]
                        sel_mol_id = sel_atom_to_mol[sel_atom_traj_id]
                    
                    # Only include if different molecule
                    if ref_mol_id != sel_mol_id:
                        filtered_indices.append(idx)
                neighbor_indices = np.array(filtered_indices)
            
            if len(neighbor_indices) > 0:
                neighbor_deltas = deltas_full[ref_atom_id, neighbor_indices, :]
                neighbor_dists = distance_matrix[ref_atom_id, neighbor_indices]
                
                neighbors = [(neighbor_deltas[j], neighbor_dists[j]) for j in range(len(neighbor_indices))]
                
                # Build molecule-level neighbor tracking for this ref atom (same as V0)
                if ref_mol_id not in mol_neighbors:
                    mol_neighbors[ref_mol_id] = {}
                
                for j in range(len(neighbor_indices)):
                    if args_dict['sel_mol'][i]:
                        # Using COM: neighbor_indices[j] is position index, map to mol_id
                        sel_mol_id = sel_posidx_to_molid[neighbor_indices[j]]
                    else:
                        # Using atoms: get molecule from atom index
                        sel_atom_traj_id = sel_indices[neighbor_indices[j]]
                        sel_mol_id = sel_atom_to_mol[sel_atom_traj_id]
                    
                    dist = neighbor_dists[j]
                    delta = neighbor_deltas[j]
                    sel_mol_id_in_pair = sorted_all_sel_mols.index(sel_mol_id) + 1
                    sel_atomID = sel_mol_first_atom[sel_mol_id] + 1
                    
                    # Track closest distance for each ref-sel molecule pair (same format as V0)
                    if sel_mol_id not in mol_neighbors[ref_mol_id]:
                        mol_neighbors[ref_mol_id][sel_mol_id] = {
                            'dist': dist,
                            'delta': delta,
                            'sel_mol_id_in_pair': sel_mol_id_in_pair,
                            'sel_atomID': sel_atomID
                        }
                    elif dist < mol_neighbors[ref_mol_id][sel_mol_id]['dist']:
                        mol_neighbors[ref_mol_id][sel_mol_id] = {
                            'dist': dist,
                            'delta': delta,
                            'sel_mol_id_in_pair': sel_mol_id_in_pair,
                            'sel_atomID': sel_atomID
                        }
                
                # VECTORIZED: Compute Qn with batch harmonics
                q6 = compute_q6_vectorized(neighbor_deltas, l=args_dict['order'][i])
            else:
                neighbors = []
                q6 = np.nan
            
            global_ref_id = len(ref_indices) * i + ref_atom_id
            if global_ref_id not in q6_per_ref_frame:
                q6_per_ref_frame[global_ref_id] = []
            q6_per_ref_frame[global_ref_id].append((q6, len(neighbors)))
            
            pair_q6_values.append(q6)
            pair_nn_values.append(len(neighbors))
        
        # Create neighbor details - one entry per reference molecule (EXACT copy of V0 logic)
        for ref_mol_id_global in sorted(mol_neighbors.keys()):
            ref_mol_id_in_pair = sorted_all_ref_mols.index(ref_mol_id_global) + 1
            ref_atomID = ref_mol_first_atom[ref_mol_id_global] + 1
            neighbors_list = []
            for sel_mol_id_global in sorted(mol_neighbors[ref_mol_id_global].keys()):
                info = mol_neighbors[ref_mol_id_global][sel_mol_id_global]
                neighbors_list.append((info['sel_mol_id_in_pair'], info['sel_atomID'], info['dist'], info['delta']))
            
            neighbor_details.append({
                'frame': frame_idx,
                'pair_idx': i,
                'ref_name': ref_name,
                'ref_id': ref_mol_id_in_pair,
                'ref_atom_id': ref_atomID,
                'sel_name': sel_name,
                'neighbors': neighbors_list
            })
        
        # Calculate per-pair average Q6 and average NN
        pair_q6_array = np.array(pair_q6_values)
        q6_by_pair[i] = np.nanmean(pair_q6_array)
        nn_by_pair[i] = np.mean(pair_nn_values) if pair_nn_values else 0.0
    
    return frame_idx, q6_by_pair, nn_by_pair, neighbor_details, q6_per_ref_frame

def read_trajectory_frames(args):
    """Read all trajectory frames."""
    frames = []
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
                    frames.append((coords, box, atom_info, frame_number))
                frame_number += 1
                for _ in range(args.skip - 1):
                    read_gro_frame(f, skip=True)
                    frame_number += 1
        else:
            for _ in range(args.end - args.begin):
                coords, box, atom_info = read_gro_frame(f, skip=False)
                if coords is None:
                    break
                if frame_number % max(1, args.skip) == 0:
                    frames.append((coords, box, atom_info, frame_number))
                frame_number += 1
    
    return frames

def calculate_q6_trajectory_combined(args, groups, is_self_pair, groups_mol, box_init, nproc=None):
    """Calculate Q6 using vectorized + parallel approach (per-pair)."""
    
    args_dict = {
        'ref': args.ref,
        'sel': args.sel,
        'ref_mol': args.ref_mol,
        'sel_mol': args.sel_mol,
        'rcut': args.rcut,
    }
    
    logging.info(f"Reading trajectory frames...")
    frames = read_trajectory_frames(args)
    logging.info(f"Read {len(frames)} frames")
    
    frame_tasks = []
    for coords, box, atom_info, frame_idx in frames:
        frame_tasks.append((coords, box, atom_info, args_dict, groups, is_self_pair, groups_mol, frame_idx))
    
    # Track per-pair time series
    npairs = len(args.ref)
    q6_timeseries_by_pair = {i: [] for i in range(npairs)}
    nn_timeseries_by_pair = {i: [] for i in range(npairs)}
    q6_per_ref = {}
    all_neighbor_details = []
    
    if nproc is None:
        nproc = os.cpu_count()
    logging.info(f"Using {nproc} processes with VECTORIZED operations")
    
    with Pool(processes=nproc) as pool:
        results = pool.map(process_frame_q6_combined, frame_tasks)
    
    # Aggregate results by pair
    for frame_idx, q6_by_pair, nn_by_pair, neighbor_details, q6_per_ref_frame in results:
        for pair_idx in range(npairs):
            q6_timeseries_by_pair[pair_idx].append((frame_idx, q6_by_pair[pair_idx]))
            nn_timeseries_by_pair[pair_idx].append((frame_idx, nn_by_pair[pair_idx]))
        all_neighbor_details.extend(neighbor_details)
        for ref_id, q6_list in q6_per_ref_frame.items():
            if ref_id not in q6_per_ref:
                q6_per_ref[ref_id] = []
            q6_per_ref[ref_id].extend(q6_list)
    
    # Sort and convert to arrays for each pair
    for pair_idx in range(npairs):
        q6_timeseries_by_pair[pair_idx].sort(key=lambda x: x[0])
        q6_timeseries_by_pair[pair_idx] = np.array([q6 for _, q6 in q6_timeseries_by_pair[pair_idx]])
        
        nn_timeseries_by_pair[pair_idx].sort(key=lambda x: x[0])
        nn_timeseries_by_pair[pair_idx] = np.array([nn for _, nn in nn_timeseries_by_pair[pair_idx]])
    
    neighbor_stats = {'n_neighbors': [], 'n_ref_atoms': [len(q6_per_ref)]}
    
    return q6_timeseries_by_pair, nn_timeseries_by_pair, q6_per_ref, neighbor_stats, all_neighbor_details

def write_q6_output(args, q6_timeseries_by_pair, nn_timeseries_by_pair, all_neighbor_details):
    """Write output - single multi-column Qn file + separate NN file."""
    npairs = len(q6_timeseries_by_pair)
    
    if npairs == 0 or len(q6_timeseries_by_pair[0]) == 0:
        logging.warning("No frames processed.")
        return
    
    base_name = os.path.splitext(args.output_file)[0]
    
    # Write Qn file with multiple columns (one per pair)
    q6_file = args.output_file
    q6_header = "# Time(ps)" + "".join([f"    Q{args.order[i]}({args.ref[i]}-{args.sel[i]})" for i in range(npairs)])
    q6_output_data = []
    for frame_idx in range(len(q6_timeseries_by_pair[0])):
        time_ps = frame_idx * args.skip
        row = [time_ps]
        for pair_idx in range(npairs):
            row.append(q6_timeseries_by_pair[pair_idx][frame_idx])
        q6_output_data.append(row)
    
    q6_array = np.array(q6_output_data)
    fmt_str = f"{'%10.3f'} {' %10.6f' * npairs}"
    np.savetxt(q6_file, q6_array, header=q6_header, fmt=fmt_str, comments='')
    logging.info(f"Qn time series written to {q6_file}")
    
    # Print statistics for each pair
    logging.info("\n--- Qn Statistics ---")
    for pair_idx in range(npairs):
        pair_name = f"{args.ref[pair_idx]}-{args.sel[pair_idx]}"
        order = args.order[pair_idx]
        q6_data = q6_timeseries_by_pair[pair_idx]
        nn_data = nn_timeseries_by_pair[pair_idx]
        
        logging.info(f"\n{pair_name} (Q{order}):")
        logging.info(f"  Mean Q{order}:     {np.nanmean(q6_data):.6f}")
        logging.info(f"  Std Dev Q{order}:  {np.nanstd(q6_data):.6f}")
        logging.info(f"  Min Q{order}:      {np.nanmin(q6_data):.6f}")
        logging.info(f"  Max Q{order}:      {np.nanmax(q6_data):.6f}")
        logging.info(f"  Mean Avg_NN: {np.nanmean(nn_data):.2f}")
    
    # Write separate NN file
    nn_file = f"{base_name}_NN.xvg"
    nn_header = "# Time(ps)" + "".join([f"    NN({args.ref[i]}-{args.sel[i]})" for i in range(npairs)])
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
    logging.info(f"Average NN time series written to {nn_file}")
    
    if args.detail:
        detail_file = os.path.splitext(args.output_file)[0] + '_detail.txt'
        with open(detail_file, 'w') as f:
            f.write("# Frame  Time(ps)  Pair_id  Ref_id  Ref_AtomID  Ref_Name  NN_id  Sel_id  Sel_AtomID  Sel_Name  Distance(nm)  Vec_X  Vec_Y  Vec_Z  Total_NN\n")
            for detail in all_neighbor_details:
                frame = detail['frame']
                time_ps = frame * args.skip
                pair_id = detail['pair_idx'] + 1
                ref_id = detail['ref_id']
                ref_atom_id = detail['ref_atom_id']
                ref_name = detail['ref_name']
                sel_name = detail['sel_name']
                total_nn = len(detail['neighbors'])
                for nn_idx, (sel_id, sel_atom_id, distance, vector) in enumerate(detail['neighbors'], start=1):
                    f.write(f"{frame:6d} {time_ps:10.3f} {pair_id:8d} {ref_id:8d} {ref_atom_id:11d} {ref_name:8s} {nn_idx:6d} {sel_id:8d} {sel_atom_id:11d} {sel_name:8s} {distance:13.6f} {vector[0]:10.6f} {vector[1]:10.6f} {vector[2]:10.6f} {total_nn:9d}\n")
        logging.info(f"Detailed information written to {detail_file}")

def main():
    """Main function."""
    start_time = time.time()
    args = get_cli_args()
    
    nref = len(args.ref)
    if args.ref_mol is None:
        args.ref_mol = [False] * nref
    elif len(args.ref_mol) == 1:
        args.ref_mol = args.ref_mol * nref
    
    nsel = len(args.sel)
    if nsel != nref:
        raise ValueError(f"Selection count must match reference count")
    if args.sel_mol is None:
        args.sel_mol = [False] * nsel
    elif len(args.sel_mol) == 1:
        args.sel_mol = args.sel_mol * nsel
    
    npair = nref
    
    # Validate and expand order argument
    if args.order is None:
        args.order = [6] * npair  # Default to Q6
    elif len(args.order) == 1:
        args.order = args.order * npair  # Expand single value to all pairs
    elif len(args.order) != npair:
        raise ValueError(f"--order must have length 1 or {npair}, got {len(args.order)}")
    
    need_mass = any(args.ref_mol) or any(args.sel_mol)
    if need_mass and args.mass_file is None:
        raise ValueError("Mass file required for COM")
    
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
    logging.info("Qn STEINHARDT - V3 COMBINED (Vectorized + Parallel)")
    logging.info("=" * 70)
    logging.info(f"Using {args.nproc or os.cpu_count()} processes with vectorized ops")
    logging.info("Command-line parameters:")
    for k, v in vars(args).items():
        if k not in ['ref', 'sel', 'ref_mol', 'sel_mol', 'rcut']:
            logging.info(f"  {k}: {v}")
    
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
    
    groups_mol = None
    if need_mass:
        logging.info("\nParsing mass file...")
        masses = parse_mass_file(args.mass_file)
        logging.info(f"Found {len(masses)} molecular groups in mass file")
    
    logging.info("\nReading structure file...")
    logging.info(f"  File: {args.structure_file}")
    with open(args.structure_file, 'r') as f:
        coords_init, box_init, atom_info_init = read_gro_frame(f)
    
    if coords_init is None:
        logging.error("Failed to read structure file.")
        sys.exit(1)
    
    logging.info(f"  Total atoms: {len(coords_init)}")
    logging.info(f"  Box dimensions: {box_init[0]:.4f} x {box_init[1]:.4f} x {box_init[2]:.4f} nm")
    
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
    
    if args.rcut is None or len(args.rcut) == 0:
        args.rcut = [np.min(box_init[:3]) / 2.0] * npair
    elif len(args.rcut) == 1:
        args.rcut = args.rcut * npair
    
    if not isinstance(args.rcut, list):
        args.rcut = list(args.rcut)
    
    logging.info("\nCutoff Distance (Rcut):")
    for i, (ref_name, sel_name) in enumerate(zip(args.ref, args.sel)):
        logging.info(f"  Pair {i+1} ({ref_name}-{sel_name}): {args.rcut[i]:.4f} nm")
    
    is_self_pair = np.zeros(npair, dtype=bool)
    logging.info("\nPair Configuration:")
    for i, (ref_name, sel_name) in enumerate(zip(args.ref, args.sel)):
        if ref_name == sel_name:
            is_self_pair[i] = True
        self_pair_str = "YES (self-comparison)" if is_self_pair[i] else "NO (cross-pair)"
        logging.info(f"  Pair {i+1}: {ref_name} <-> {sel_name}")
        logging.info(f"    Self-pair: {self_pair_str}")
        logging.info(f"    COM mode: ref={args.ref_mol[i]}, sel={args.sel_mol[i]}")
    
    logging.info("\n--- Calculation ---")
    logging.info("Processing trajectory... [VECTORIZED + PARALLEL]")
    q6_timeseries_by_pair, nn_timeseries_by_pair, q6_per_ref, neighbor_stats, all_neighbor_details = calculate_q6_trajectory_combined(
        args, groups, is_self_pair, groups_mol, box_init, nproc=args.nproc
    )
    
    logging.info("\n--- Output ---")
    write_q6_output(args, q6_timeseries_by_pair, nn_timeseries_by_pair, all_neighbor_details)
    
    logging.info("\n" + "=" * 70)
    logging.info("CALCULATION COMPLETED")
    logging.info("=" * 70)
    
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
    logging.info(f"Processes:     {args.nproc or os.cpu_count()}")

if __name__ == '__main__':
    main()
