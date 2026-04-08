# Qn Steinhardt Order Parameter Calculator

A generalized implementation of the Steinhardt order parameter for analyzing local crystalline order in molecular dynamics simulations using GROMACS trajectories.

## Overview

This project provides two versions of the Qn (Steinhardt) order parameter calculator:
- **V0**: Sequential processing - suitable for single-pair analysis or testing
- **V3**: Vectorized + Parallelized - optimized for multi-pair analysis and large trajectories

### Key Features

- **Flexible Order Parameter**: Specify different order parameters (Q1, Q2, Q4, Q6, Q11, etc.) for different pairs
- **Multi-pair Analysis**: Compute order parameters for multiple reference-neighbor pairs simultaneously
- **Center-of-Mass (COM) Support**: Optional COM calculation for molecular systems
- **Adaptive Cutoff**: Automatic or manual definition of neighbor cutoff distances
- **Detailed Neighbor Logging**: Track individual neighbor atoms/molecules with distances
- **Performance Optimized**: V3 uses vectorized NumPy operations and parallel processing

## Installation

### Requirements

- Python 3.7+
- NumPy
- SciPy (Bessel functions)

### Setup

```bash
pip install numpy scipy
```

## Usage

### Basic Command Structure

```bash
python3 Qn_order_parameter.py \
  -s structure.gro \
  -f trajectory.gro \
  -n index.ndx \
  --ref GROUP1 [GROUP2 ...] \
  --sel GROUP1 [GROUP2 ...] \
  [--order O1 O2 ...] \
  [--rcut CUTOFF] \
  [OPTIONS]
```

### V0 vs V3 Selection

**Use V0 (Sequential) when:**
- Analyzing single or few pairs
- Debugging or validating results
- Memory is limited

**Use V3 (Parallel) when:**
- Analyzing multiple pairs (3+)
- Processing long trajectories
- Performance is critical

```bash
# V0: Sequential version
python3 Qn_order_parameter.py ...

# V3: Parallel version (4 processes)
python3 Qn_order_parameter_v3_combined.py -j 4 ...
```

## Input Files Required

1. **Structure file** (`*.gro`): GROMACS structure format defining the system
2. **Trajectory file** (`*.gro`): GROMACS trajectory in .gro format
3. **Index file** (`*.ndx`): GROMACS index file defining atom groups
4. **Mass file** (optional): Atomic masses for COM calculation (if `--ref_mol` or `--sel_mol` used)

### Index File Format (*.ndx)

```
[ System ]
    1    2    3    4 ...

[ SOL ]
    5    6    7    8    9   10

[ ION ]
   11   12
```

### Mass File Format

Plain text with residue/atom names and masses:
```
SOL   18.015
ION   39.098
```

## Output Files

| File | Description |
|------|-------------|
| `Qn.xvg` | Time series of Qn values for each pair |
| `Qn_NN.xvg` | Average number of neighbors over time |
| `Qn.log` | Detailed execution log with statistics |
| `Qn_detail.txt` | (if `--detail`) Neighbor list for each molecule |

## Complete Examples

### Example 1: Single Pair, Default Q6

```bash
python3 Qn_order_parameter.py \
  -s simulation.gro \
  -f trajectory.gro \
  -n index.ndx \
  --ref WATER \
  --sel IONS
```

**Input:**
- Structure: `simulation.gro` (GROMACS structure file)
- Trajectory: `trajectory.gro` (100 frames)
- Index file with groups: `WATER` (500 atoms) and `IONS` (50 atoms)

**Output:**
```
Time(ps)      Q6(WATER-IONS)
    0.0           0.521845
    1.0           0.523127
    2.0           0.519834
    ...
```

### Example 2: Multiple Pairs with Different Orders

```bash
python3 Qn_order_parameter.py \
  -s system.gro \
  -f trajectory.gro \
  -n index.ndx \
  --ref SOL ION POLYMER \
  --sel SOL ION POLYMER \
  --order 6 4 2 \
  --rcut 0.5 0.6 0.7
```

**Output Header:**
```
# Time(ps)       Q6(SOL-SOL)      Q4(ION-ION)    Q2(POLYMER-POLYMER)
      0.0          0.521845          0.412134            0.187204
      1.0          0.523127          0.415892            0.189467
```

### Example 3: Center-of-Mass Calculation

```bash
python3 Qn_order_parameter.py \
  -s molecules.gro \
  -f trajectory.gro \
  -n index.ndx \
  --ref WATER \
  --sel IONS \
  --ref_mol true \
  --sel_mol false \
  --mass_file masses.dat \
  --rcut 1.0
```

- `--ref_mol true`: Use COM of WATER molecules
- `--sel_mol false`: Use individual ION atoms
- `masses.dat` contains atomic masses
- Cutoff: 1.0 nm for all pairs

### Example 4: Parallel Processing with Detailed Output

```bash
python3 Qn_order_parameter_v3_combined.py \
  -s system.gro \
  -f trajectory.gro \
  -n index.ndx \
  --ref GROUP1 GROUP2 GROUP3 \
  --sel GROUP1 GROUP2 GROUP3 \
  --order 6 6 6 \
  -j 4 \
  --detail \
  --begin 0 \
  --end 100 \
  --skip 2
```

**Options:**
- `-j 4`: Use 4 parallel processes
- `--detail`: Write detailed neighbor information
- `--begin 0`: Start from frame 0
- `--end 100`: Process up to frame 100
- `--skip 2`: Process every 2nd frame

### Example 5: Single Order for All Pairs

```bash
python3 Qn_order_parameter.py \
  -s sim.gro \
  -f traj.gro \
  -n index.ndx \
  --ref A B C D \
  --sel X Y Z W \
  --order 6
```

Single value automatically expands: `--order 6` → `[6, 6, 6, 6]` for 4 pairs

## Command-line Arguments

### Required Arguments

| Argument | Description |
|----------|-------------|
| `-s`, `--structure_file` | Structure file (.gro) |
| `-f`, `--traj_file` | Trajectory file (.gro) |
| `-n`, `--index_file` | Index file (.ndx) |
| `--ref` | Reference group name(s) |
| `--sel` | Selection group name(s) |

### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `-o`, `--output_file` | `Qn.xvg` | Output filename |
| `--order` | `6` for all pairs | Qn order parameter(s) |
| `--rcut` | `min(box)/2` | Cutoff distance (nm) |
| `--ref_mol` | `false` | Use COM for reference |
| `--sel_mol` | `false` | Use COM for selection |
| `--mass_file` | None | Mass file for COM |
| `--begin` | `0` | First frame to process |
| `--end` | `-1` (last) | Last frame to process |
| `--skip` | `1` | Process every Nth frame |
| `--detail` | False | Write detailed neighbor output |
| `--debug` | False | Enable debug logging |
| `-j`, `--nproc` | (V3 only) | Number of processes |

## Output Explanation

### Qn.xvg (Main Output)

Multi-column time series file:
```
# Time(ps)       Q6(REF1-SEL1)      Q4(REF2-SEL2)      Q2(REF3-SEL3)
      0.000           0.520184            0.412345            0.195623
      1.000           0.523847            0.415892            0.198134
      2.000           0.521956            0.410123            0.192847
```

### Qn_NN.xvg (Neighbor Counts)

Average number of neighbors per reference atom/molecule:
```
# Time(ps)          NN(REF1-SEL1)       NN(REF2-SEL2)       NN(REF3-SEL3)
      0.000              12.5                8.3                 15.2
      1.000              12.7                8.1                 15.4
      2.000              12.4                8.5                 15.1
```

### Qn.log (Execution Log)

Contains:
- System parameters
- Processing statistics
- Per-pair statistics (mean, std dev, min, max)
- Execution time

```
======================================================================
Qn STEINHARDT ORDER PARAMETER CALCULATION
======================================================================
Command-line parameters:
  structure_file: simulation.gro
  traj_file: trajectory.gro
  index_file: index.ndx
  ...

Group Summary:
  Pair 1: WATER <-> IONS
    Atoms in WATER: 1500 (3 atoms/molecule = 500 molecules)
    Atoms in IONS: 100
    Self-pair: False
    Cutoff: 0.5000 nm

--- Qn Statistics ---

WATER-IONS (Q6):
  Mean Q6:     0.521456
  Std Dev Q6:  0.012345
  Min Q6:      0.489123
  Max Q6:      0.556789
  Mean Avg_NN: 12.34
```

### Qn_detail.txt (Optional, with --detail)

Detailed neighbor information for each reference molecule/atom:

```
Frame 0
-------
Pair 0: WATER <-> IONS
  Molecule 1 (atom 1):
    Total_NN: 13
    Distances: [0.234, 0.287, 0.341, ...]
    Neighbors: [(atom_id=45, sel_mol_id=2, dist=0.234), ...]
  
  Molecule 2 (atom 4):
    Total_NN: 11
    Distances: [0.256, 0.312, ...]
    ...
```

## Scientific Background

The Steinhardt order parameter Qn characterizes local crystalline ordering:

$$Q_n = \sqrt{\frac{4\pi}{2n+1} \sum_{m=-n}^{n} |\bar{q}_n^m|^2}$$

where:
- $q_n^m(\mathbf{r}_i) = Y_n^m(\theta, \phi)$ is the spherical harmonic
- $\bar{q}_n^m$ is the average over neighboring atoms
- Common values: Q6 for hexatic order, Q4 for tetrahedral, Q2 for orientational

### Reference

Steinhardt, P. J., Nelson, D. R., & Ronchetti, M. (1981).
"Bond-orientational order in liquids and glasses."
*Physical Review B*, 28(2), 784–805.

## Performance Notes

### V0 (Sequential)
- Speed: ~1-2 frames/second per pair
- Memory: ~100-200 MB
- Best for: Development, debugging, single pairs

### V3 (Parallel + Vectorized)
- Speed: ~5-15 frames/second per pair (3-4x from vectorization, 3-4x from parallelization)
- Memory: ~200-400 MB (scales with process count)
- Best for: Production runs, multiple pairs, long trajectories

### Optimization Tips

1. **Cutoff distance**: Larger cutoffs increase computation but improve statistics
2. **Frame skipping**: Use `--skip 2` or `--skip 5` for long trajectories
3. **Parallel processes**: Optimal is usually CPU count or CPU count - 1
4. **Frame range**: Process only frames of interest with `--begin` and `--end`

## Troubleshooting

### Common Issues

**Error: "Index group not found"**
- Check group name spelling in index file (case-sensitive)
- Verify index file format is correct

**Error: "Mass file required but not provided"**
- Add `--mass_file` when using `--ref_mol` or `--sel_mol`
- Verify mass file has correct format

**Slow performance**
- Try V3 version with parallel processing
- Increase `--skip` value for long trajectories
- Reduce frame range with `--begin` and `--end`

**NaN values in output**
- Usually indicates no neighbors within cutoff
- Try increasing `--rcut`
- Verify selection groups have atoms

## File Structure

```
├── Qn_order_parameter.py          # V0: Sequential version
├── Qn_order_parameter_v3_combined.py  # V3: Parallel version
├── README.md                       # This file
├── MANUAL.md                       # Technical documentation
└── examples/
    ├── example1_single_pair/
    ├── example2_multi_pair/
    └── example3_with_com/
```

## License

[Your License Here]

## Contact

For questions or issues, contact: [Your Email/GitHub]
