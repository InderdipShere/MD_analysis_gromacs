# Quick Start Guide

## 1. Installation

```bash
# Clone or download the repository
cd your_repository

# Install dependencies
pip install numpy scipy
```

## 2. Prepare Your Input Files

Create three essential files:
- `structure.gro` - GROMACS structure
- `trajectory.gro` - GROMACS trajectory
- `index.ndx` - GROMACS index file

See `examples/README_INPUT_FILES.md` for details.

## 3. Run Basic Analysis

```bash
# Single pair with default Q6
python3 Qn_order_parameter.py \
  -s structure.gro \
  -f trajectory.gro \
  -n index.ndx \
  --ref WATER \
  --sel IONS
```

## 4. Check Results

```bash
# View output
cat Qn.xvg

# Check statistics in log
cat Qn.log
```

## Files in This Repository

### Main Scripts
- `Qn_order_parameter.py` - Sequential version (V0)
- `Qn_order_parameter_v3_combined.py` - Parallel version (V3)

### Documentation
- `README.md` - Full user guide with examples
- `MANUAL.md` - Technical details and algorithm explanation
- `QUICKSTART.md` - This file
- `examples/EXAMPLES.md` - Input/output examples

### Examples
- `examples/sample_structure.gro` - Sample structure file
- `examples/sample_index.ndx` - Sample index file
- `examples/sample_masses.dat` - Sample mass file
- `examples/README_INPUT_FILES.md` - Input file format guide

## Common Commands

### Default Analysis (Q6 for all pairs)

```bash
python3 Qn_order_parameter.py \
  -s structure.gro -f trajectory.gro -n index.ndx \
  --ref GROUP1 --sel GROUP2
```

### Multiple Pairs with Different Orders

```bash
python3 Qn_order_parameter.py \
  -s structure.gro -f trajectory.gro -n index.ndx \
  --ref G1 G2 G3 --sel S1 S2 S3 \
  --order 6 4 2
```

### Parallel Processing (Faster)

```bash
python3 Qn_order_parameter_v3_combined.py \
  -s structure.gro -f trajectory.gro -n index.ndx \
  --ref GROUP1 --sel GROUP2 \
  -j 4
```

### With Center-of-Mass

```bash
python3 Qn_order_parameter.py \
  -s structure.gro -f trajectory.gro -n index.ndx \
  --ref GROUP1 --sel GROUP2 \
  --ref_mol true --sel_mol false \
  --mass_file masses.dat
```

## Output Files

| File | Description |
|------|-------------|
| `Qn.xvg` | Time series of order parameters |
| `Qn_NN.xvg` | Average neighbor counts |
| `Qn.log` | Detailed analysis log |
| `Qn_detail.txt` | (optional) Neighbor lists |

## Troubleshooting

**Error: "Group not found"**
```bash
grep "\[" index.ndx  # Check group names
```

**Error: "Module numpy not found"**
```bash
pip install numpy scipy
```

**Error: "All NaN values"**
```bash
# Try larger cutoff
--rcut 1.0
```

## Next Steps

1. Read `README.md` for complete documentation
2. Check `examples/EXAMPLES.md` for detailed examples
3. Review `MANUAL.md` for technical details
4. Run your own analysis!

## Need Help?

See:
- `README.md` - Comprehensive guide
- `MANUAL.md` - Technical implementation details
- `examples/` - Example commands and files
- Script help: `python3 Qn_order_parameter.py --help`
