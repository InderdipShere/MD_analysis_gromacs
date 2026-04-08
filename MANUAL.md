# Technical Manual - Qn Steinhardt Order Parameter

## Table of Contents

1. [Algorithm Overview](#algorithm-overview)
2. [Implementation Details](#implementation-details)
3. [Code Architecture](#code-architecture)
4. [Key Functions](#key-functions)
5. [Performance Considerations](#performance-considerations)
6. [Molecular vs Atomic Modes](#molecular-vs-atomic-modes)
7. [Output File Formats](#output-file-formats)

## Algorithm Overview

### Theoretical Background

The Steinhardt order parameter Qn characterizes local crystalline structure by measuring the correlation of bond angles around each particle:

#### Formula

$$Q_n = \sqrt{\frac{4\pi}{2n+1} \sum_{m=-n}^{n} |\bar{q}_n^m|^2}$$

where:

$$q_n^m(i) = \frac{1}{N_b} \sum_{j=1}^{N_b} Y_n^m(\theta_{ij}, \phi_{ij})$$

### Computational Steps

1. **Neighbor Detection**: For each reference atom/molecule, find all selection atoms/molecules within cutoff distance
2. **Spherical Harmonic Calculation**: Compute $Y_n^m(\theta, \phi)$ for each neighbor pair
3. **Harmonic Averaging**: Average harmonics over all neighbors
4. **Order Parameter**: Calculate magnitude of averaged harmonics

## Implementation Details

### Core Components

#### 1. Spherical Harmonics Computation

```python
def compute_spherical_harmonics(theta, phi, l=6):
    """
    Compute all spherical harmonics Y_l^m(theta, phi) for m = -l to l.
    
    Parameters:
    -----------
    theta : array-like
        Polar angle (0 to π)
    phi : array-like
        Azimuthal angle (0 to 2π)
    l : int
        Harmonic order
    
    Returns:
    --------
    complex_array of shape (..., 2*l+1)
        Harmonics for m = -l to l
    """
```

**Implementation:**
- Uses `scipy.special.sph_harm_y()` for accurate computation
- Handles vectorized input for efficiency (V3)
- Converts Cartesian deltas (dx, dy, dz) to spherical coordinates

#### 2. Neighbor Filtering

**Algorithms:**
1. Distance calculation using `scipy.spatial.distance.cdist()`
2. Apply cutoff threshold
3. **Molecule-level self-pair exclusion**: If ref and sel are same group, exclude neighbors in same molecule

**Key Logic:**
```python
# Molecular mode: map atom indices to molecule IDs
ref_atom_to_mol = build_molecule_mapping(atom_info, ref_indices)

# For self-pairs, exclude atoms in same molecule
if is_self_pair[pair_idx]:
    filtered_neighbors = [
        atom for atom in neighbors 
        if ref_mol_id != sel_atom_to_mol[atom]
    ]
```

#### 3. Molecule Mapping

**Residue-based:** Each unique residue number represents one molecule
- Sorted residue numbers: ensures consistent numbering across frames
- Molecule ID = position in sorted residue list
- First atom of each residue stored for reference

```python
unique_residues = sorted(set(atom_info[idx][0] for idx in group_indices))
residue_to_mol_id = {res: i+1 for i, res in enumerate(unique_residues)}
```

### Position Modes

#### Atomic Mode (`--ref_mol false`, `--sel_mol false`)

Each atom is treated individually:
- Reference position: individual atom position
- Selection neighbors: individual atoms within cutoff
- Numbering: atoms 1, 2, 3...

#### Center-of-Mass Mode (`--ref_mol true` or `--sel_mol true`)

Molecules treated as single points:
- Reference position: COM of reference molecule
- Selection neighbors: COM of selection molecules
- Numbering: molecules 1, 2, 3... (one per residue)

**COM Calculation:**
```python
com = sum(mass[i] * position[i]) / sum(mass[i])
```

## Code Architecture

### V0: Sequential Version

```
main()
  ├─ parse_cli_args()
  ├─ read_structure()
  ├─ parse_ndx()
  ├─ read_trajectory()
  ├─ validate_arguments()
  └─ for each frame:
      └─ process_frame_q6()
          ├─ build_molecule_mapping()
          ├─ compute distance matrix
          ├─ for each reference atom/mol:
          │   ├─ find neighbors
          │   ├─ filter by cutoff
          │   ├─ exclude self-pairs
          │   └─ compute_Q6_for_reference()
          │       ├─ compute_spherical_harmonics()
          │       ├─ average over neighbors
          │       └─ return Qn value
          └─ collect results
  └─ write_q6_output()
```

### V3: Parallel + Vectorized Version

```
main()
  ├─ parse_cli_args()
  ├─ read_structure()
  ├─ prepare_frames_for_pool()
  └─ Pool.map(process_frame_q6_combined, frame_list)
      └─ worker: process_frame_q6_combined()
          ├─ VECTORIZED: cdist() for all distances
          ├─ for each pair:
          │   ├─ vectorized filtering
          │   ├─ VECTORIZED: compute_q6_vectorized()
          │   │   ├─ batch spherical harmonic calc
          │   │   ├─ vectorized averaging
          │   │   └─ return Qn array
          │   └─ process results by molecule
          └─ return frame results
  └─ collect from all workers + write_q6_output()
```

**Parallelization Strategy:**
- Frame-level parallelization: each worker processes one frame
- All pairs computed for a frame in same worker → better cache locality
- Results collected and written sequentially

## Key Functions

### compute_Q6_for_reference (V0)

```python
def compute_Q6_for_reference(neighbors_data, l=6):
    """
    Compute order parameter for single reference.
    
    Parameters:
    -----------
    neighbors_data : list of (delta_vector, distance)
        Displacement vectors and distances to neighbors
    l : int
        Harmonic order
    
    Returns:
    --------
    float : Qn value in range [0, 1]
    """
    if len(neighbors_data) == 0:
        return np.nan
    
    # Convert deltas to spherical coordinates
    deltas_arr = np.array([d[0] for d in neighbors_data])
    theta, phi = convert_to_spherical(deltas_arr)
    
    # Compute harmonics for all neighbors
    harmonics = compute_spherical_harmonics(theta, phi, l)  # shape: (N_neighbors, 2*l+1)
    
    # Average over neighbors
    bar_qlm = np.mean(harmonics, axis=0)
    
    # Compute magnitude
    ql = np.sqrt(4*np.pi/(2*l+1) * np.sum(np.abs(bar_qlm)**2))
    
    return ql
```

### build_molecule_mapping

```python
def build_molecule_mapping(atom_info, indices):
    """
    Create atom-to-molecule mapping based on residue numbers.
    
    Parameters:
    -----------
    atom_info : ndarray
        Shape (N_atoms, 6): [res_nr, res_name, atom_name, atom_nr, x, y]
    indices : list
        Atom indices in this group
    
    Returns:
    --------
    atom_to_mol : dict
        Maps trajectory atom ID → molecule ID (1-indexed)
    mol_first_atom : dict
        Maps molecule ID → first atom index in trajectory
    """
    # Extract unique residues (sorted for consistency)
    unique_residues = sorted(set(int(atom_info[idx][0]) for idx in indices))
    
    # Create mapping: residue_number → 1, 2, 3, ...
    residue_to_mol_id = {res: i+1 for i, res in enumerate(unique_residues)}
    
    # Map atoms
    atom_to_mol = {}
    mol_first_atom = {}
    for traj_idx in indices:
        res_nr = int(atom_info[traj_idx][0])
        mol_id = residue_to_mol_id[res_nr]
        atom_to_mol[traj_idx] = mol_id
        
        if mol_id not in mol_first_atom:
            mol_first_atom[mol_id] = traj_idx
    
    return atom_to_mol, mol_first_atom
```

### process_frame_q6 (V0)

Processes single frame:
1. Build molecule mappings for each pair
2. Compute distances (cdist for all pairs at once)
3. For each pair:
   - For each reference mol/atom:
     - Get neighbors within cutoff
     - Filter self-pairs
     - Compute Qn
4. Return per-pair results

```python
def process_frame_q6(args, groups, is_self_pair, groups_mol, coords, box, 
                     atom_info, q6_per_ref, frame_idx):
    """Process single frame using sequential algorithm."""
    
    npairs = len(args.ref)
    q6_by_pair = {}
    nn_by_pair = {}
    
    for pair_idx, (ref_name, sel_name) in enumerate(zip(args.ref, args.sel)):
        ref_indices = groups[ref_name]
        sel_indices = groups[sel_name]
        
        # Setup molecule mappings
        ref_atom_to_mol, ref_mol_first_atom = build_molecule_mapping(atom_info, ref_indices)
        sel_atom_to_mol, sel_mol_first_atom = build_molecule_mapping(atom_info, sel_indices)
        
        # Get positions (atomic or COM)
        if args.ref_mol[pair_idx]:
            ref_positions = compute_com(coords, ref_indices, atom_info, masses)
        else:
            ref_positions = coords[ref_indices]
        
        # Compute all distances at once
        distances = cdist(ref_positions, sel_positions)
        
        # Process each reference
        for ref_idx, distances_to_all_sel in enumerate(distances):
            neighbor_indices = np.where(distances_to_all_sel < rcut)[0]
            
            # Filter self-pair if needed
            if is_self_pair[pair_idx]:
                neighbor_indices = filter_self_molecules(neighbor_indices, ...)
            
            # Compute Qn
            q6 = compute_Q6_for_reference(neighbor_deltas, l=args.order[pair_idx])
```

### process_frame_q6_combined (V3)

Vectorized version:
- Computes all distances once with cdist
- Vectorized spherical harmonic computation
- Batch processing of multiple neighbors

```python
def compute_q6_vectorized(deltas, l=6):
    """
    Compute Q values in vectorized form.
    
    Parameters:
    -----------
    deltas : ndarray
        Shape (N_neighbors, 3) - displacement vectors
    l : int
        Harmonic order
    
    Returns:
    --------
    float : Qn value
    """
    # Vectorized spherical coordinate conversion
    distances = np.linalg.norm(deltas, axis=1)
    cos_theta = deltas[:, 2] / distances
    phi = np.arctan2(deltas[:, 1], deltas[:, 0])
    theta = np.arccos(cos_theta)
    
    # Vectorized harmonic computation
    m_values = np.arange(-l, l+1)
    harmonics = np.zeros((len(deltas), 2*l+1), dtype=complex)
    for i, m in enumerate(m_values):
        harmonics[:, i] = sph_harm_y(m, l, theta, phi)
    
    # Average and magnitude
    bar_qlm = np.mean(harmonics, axis=0)
    ql = np.sqrt(4*np.pi/(2*l+1) * np.sum(np.abs(bar_qlm)**2))
    
    return ql
```

## Molecular vs Atomic Modes

### Atomic Mode (`--ref_mol false --sel_mol false`)

**Use case:** Tracking order around individual particles

```
Atom 1 -- Qn value (based on neighbors within cutoff)
Atom 2 -- Qn value
...
Atom N -- Qn value
```

**Molecule mapping still used for:**
- Self-pair exclusion (atoms from same molecule excluded)
- Neighbor deduplication (keeps only closest neighbor per molecule)

### Molecular Mode (`--ref_mol true --sel_mol true`)

**Use case:** Analyzing molecular organization

```
Molecule 1 (COM) -- Qn value (based on neighboring molecules' COMs)
Molecule 2 (COM) -- Qn value
...
Molecule M -- Qn value
```

**Features:**
- Distances computed between molecule centers
- All atoms of molecule share same Qn value
- Requires mass file for COM calculation

### Mixed Mode (`--ref_mol true --sel_mol false`)

**Example:** Tracking water coordination around ions
- Reference: Water molecules (COM-based)
- Selection: Individual ions (atomic)

## Output Analysis

### Qn Values Interpretation

| Range | Interpretation |
|-------|-----------------|
| 0.00 - 0.20 | No local order |
| 0.20 - 0.40 | Weak order |
| 0.40 - 0.60 | Moderate order |
| 0.60 - 0.80 | Strong order |
| 0.80 - 1.00 | Perfect crystalline order |

### Statistics in Log File

```
Mean Q6:     0.521456    ← Average over time and atoms
Std Dev Q6:  0.012345    ← Fluctuation magnitude
Min Q6:      0.489123    ← Most disordered frame
Max Q6:      0.556789    ← Most ordered frame
Mean Avg_NN: 12.34       ← Average coordination number
```

## Neighbor Deduplication

**Problem:** At molecular level, multiple atoms can be within cutoff

**Solution:** Keep only closest atom per candidate molecule

```
Molecule A has atoms 5, 6, 7 within cutoff
Distances: 0.35, 0.38, 0.41 nm
Decision: Keep only atom 5 (distance 0.35)

Result: One neighbor entry per unique molecule
```

This prevents:
- Artificial inflation of neighbor count
- Ballistic artifacts from molecular geometry

## Performance Optimization

### V0 Timing Breakdown

Typical 4-pair analysis (500 atoms per ref, 5000 atoms per sel):
- Distance computation: 60%
- Spherical harmonics: 30%
- I/O and overhead: 10%

### V3 Speedups

1. **Vectorization (2-3x faster)**
   - Single `cdist()` call for all pairs
   - Batch spherical harmonic computation
   - Full NumPy operations (C-accelerated)

2. **Parallelization (3-4x faster)**
   - Each CPU core processes one frame
   - Overhead: ~5-10% from multiprocessing
   - Scales well up to (CPU count - 1) workers

**Combined:** 6-12x faster for multi-pair analysis

### Memory Usage

**V0:**
- Distance matrix: `N_ref × N_sel × 8 bytes`
- Harmonics: `N_neighbors × (2*l+1) × 16 bytes`
- Total: ~100-200 MB

**V3 per worker:**
- Same as V0, but only one frame in memory
- Multiprocessing overhead: ~50-100 MB per worker
- Total: (workers × 150-300) MB

## Cutoff Distance Selection

### Strategies

1. **Auto (default):** `min(box_dimension) / 2`
   - Ensures no periodic image double-counting
   - Conservative choice

2. **Fixed:** Specify single value
   - Applied to all pairs
   - Use when physical basis known

3. **Per-pair:** Specify N values for N pairs
   - Tailored to different scales
   - More expensive computationally

### Determining Optimal Cutoff

1. Analyze radial distribution function (RDF)
2. Choose cutoff at first minimum after first peak
3. Verify results converge with larger cutoffs

## Validation & Debugging

### Enable Debug Output

```bash
python3 Qn_order_parameter.py ... --debug
```

Produces:
- Detailed coordinate conversions
- All harmonic values
- Neighbor lists per frame
- Distance matrix statistics

### Check Output Consistency

1. **V0 vs V3:** Should give identical (within floating-point precision) results
2. **Statistics:** Compare with literature values for known systems
3. **Trends:** Check for anomalies (sudden jumps, NaN values)

### Common Validity Checks

```python
# No NaN values
assert not np.isnan(qn_values).any()

# Values in valid range [0, 1]
assert (qn_values >= 0).all() and (qn_values <= 1).all()

# Statistics make sense
assert np.std(qn_values) < np.mean(qn_values)  # Usually

# Neighbor counts reasonable
assert mean_nn > 0 and mean_nn < total_atoms
```

## References

1. **Original Work:**
   - Steinhardt, P. J., Nelson, D. R., & Ronchetti, M. (1981)
   - "Bond-orientational order in liquids and glasses"
   - Physical Review B, 28(2), 784–805

2. **Extended Applications:**
   - Lechner, W., & Dellago, C. (2008)
   - "Accurate determination of crystal structures based on averaged local bond order parameters"
   - The Journal of Chemical Physics, 129(11), 114707

3. **GROMACS Tools:**
   - Official documentation: https://manual.gromacs.org/
   - Index file format: https://manual.gromacs.org/documentation/current/reference/file-formats.html#ndx

## Troubleshooting Guide

### Issue: All NaN values in output

**Causes:**
- Cutoff too small (no neighbors found)
- Index group doesn't exist
- Trajectory frames don't match structure

**Fixes:**
```bash
# Increase cutoff
--rcut 1.0

# Check groups
gmx make_ndx -f structure.gro -o index.ndx

# Verify frame compatibility
gmx check -f trajectory.gro -s structure.gro
```

### Issue: Slow performance (V0)

**Solutions:**
1. Switch to V3: ~5-10x faster for 3+ pairs
2. Reduce frame range: `--begin 0 --end 100`
3. Increase skip: `--skip 5` (process every 5th frame)
4. Reduce cutoff: Fewer neighbors = less computation

### Issue: Memory errors (V3 with many workers)

**Solutions:**
```bash
# Reduce process count
-j 2  # Instead of -j 8

# Or use V0 for large systems
python3 Qn_order_parameter.py ...
```

## Future Enhancements

Potential improvements:
1. **Bonding detection:** Use bond topology instead of cutoff
2. **Weighted harmonics:** Account for distance in weighting
3. **Symmetry parameters:** Wl (local symmetry characterization)
4. **GPU acceleration:** CUDA/OpenCL for massive speedup
