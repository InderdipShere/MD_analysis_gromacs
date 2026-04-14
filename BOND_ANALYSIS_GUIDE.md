# GROMACS Bond Analysis Tool (gmx_bond.py)

## Overview

The `gmx_bond.py` script analyzes bond formation between atoms or molecules in GROMACS molecular dynamics trajectories. It detects bonds based on distance cutoffs and provides detailed statistics including:

- **Total number of bonds** for each bond type (number of pairs)
- **Bond length statistics** (optional)
- **Support for atom-based and center-of-mass (COM) reference positions**
- **Periodic boundary condition handling** using minimum image convention
- **Vectorized calculations** for high performance
- **Detailed logging and statistics**

## Key Features

### 1. **Flexible Input**
   - Support multiple reference and selection groups
   - Explicit pairing: reference group i pairs with selection group i
   - Atom-based or center-of-mass positions

### 2. **Robust Physics**
   - Minimum image convention (MIC) for periodic boundaries
   - Proper handling of molecules split across box boundaries
   - Distance-based bond detection

### 3. **Output Format**
   - Total bond counts per frame (not normalized)
   - Statistics summary in log file
   - Optional bond length output
   - Comment lines showing total number of reference and selection atoms/molecules

### 4. **Vectorized Computation**
   - NumPy broadcasting for efficient pairwise distance calculations
   - Handles thousands of atomic pairs efficiently

---

## Installation

### Requirements
- Python 3.7+
- NumPy
- GROMACS (for trajectory files)

### Setup
```bash
pip install numpy
```

### Make executable
```bash
chmod +x gmx_bond.py
```

---

## Usage

### Basic Command Structure

```bash
python3 gmx_bond.py \
  -s structure.gro \
  -f trajectory.gro \
  -n index.ndx \
  --ref GROUP1 [GROUP2 ...] \
  --sel GROUP1 [GROUP2 ...] \
  --rcut CUTOFF1 [CUTOFF2 ...] \
  [OPTIONS]
```

### Required Arguments

| Argument | Description | Example |
|----------|-------------|---------|
| `-s, --structure_file` | Initial structure file (.gro) | `system.gro` |
| `-f, --traj_file` | Trajectory file (.gro format) | `trajectory.gro` |
| `-n, --index_file` | GROMACS index file (.ndx) | `index.ndx` |
| `--ref` | Reference group name(s) | `LI` or `LI Cl` |
| `--sel` | Selection group name(s) | `Cl` or `Cl LI` |
| `--rcut` | Cutoff distance(s) in nm | `0.35` or `0.35 0.40` |

**Important**: Number of `--ref` groups must equal number of `--sel` groups (explicit pairing).

### Optional Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `-o, --output_file` | Output file name | `bonds.xvg` |
| `--BL` | Calculate bond length statistics | False |
| `--ref_mol` | Use COM for reference groups (true/false) | False |
| `--sel_mol` | Use COM for selection groups (true/false) | False |
| `--mass_file` | File with atomic masses for COM calculation | None |
| `--begin` | First frame to analyze | 0 |
| `--end` | Last frame to analyze (-1 = end of trajectory) | -1 |
| `--skip` | Read every Nth frame | 1 |
| `--debug` | Enable debug mode | False |

---

## Input File Preparation

### 1. GROMACS Index File (.ndx)

Create or modify an index file with your desired groups using:

```bash
gmx make_ndx -f system.gro -o index.ndx
```

Example index file structure:
```
[ LI ]
  1    2    3    4    5 ...

[ Cl ]
 10   11   12   13   14 ...

[ SYSTEM ]
  1    2    3    4    5  10   11   12   13   14 ...
```

### 2. Structure File (.gro)

Either your initial structure or any frame from the trajectory:
```bash
# Used for initialization (atom names, molecule numbers)
```

### 3. Trajectory File (.gro)

GROMACS trajectory in .gro format:
```bash
# Convert from .xtc/.trr if needed
gmx trjconv -f trajectory.xtc -s system.tpr -o trajectory.gro
```

---

## Output Files

### Main Output: `bonds.xvg`

Tab-separated values with columns:
- **Column 1**: Frame number
- **Column 2+**: Total bond count for each bond type

**Example output:**
```
# frame Total_bonds_LI-Cl
0      100.000000
1      102.000000
2       98.000000
...
# Total number of Ref:
# LI: 50
# Total number of sel:
# Cl: 75
```

### Log File: `bonds.log`

Comprehensive statistics including:
- System information (hostname, platform, Python version)
- Input parameters
- Group counts
- Frame-by-frame statistics
- Mean, min, max bond counts
- Execution time

### Optional Bond Length Output: `bonds_BL.xvg`

(Only if `--BL` flag is used)
- Column 2+: Average bond length for each pair

---

## Examples

### Example 1: Simple Atom-based Bond Analysis

Analyze Li-Cl bonds in a single ionic liquid simulation:

```bash
python3 gmx_bond.py \
  -s system.gro \
  -f trajectory.gro \
  -n index.ndx \
  --ref LI \
  --sel Cl \
  --rcut 0.35 \
  -o li_cl_bonds.xvg
```

**Output interpretation:**
- Total number of Li-Cl bonds per frame
- No normalization (raw counts)

### Example 2: Multiple Bond Types

Analyze multiple bond types simultaneously:

```bash
python3 gmx_bond.py \
  -s system.gro \
  -f trajectory.gro \
  -n index.ndx \
  --ref LI Cl Na \
  --sel Cl Na LI \
  --rcut 0.35 0.40 0.38 \
  -o multi_bonds.xvg
```

**Pairs analyzed:**
- LI-Cl with cutoff 0.35 nm
- Cl-Na with cutoff 0.40 nm
- Na-LI with cutoff 0.38 nm

### Example 3: With Bond Length Statistics

Include average bond length in output:

```bash
python3 gmx_bond.py \
  -s system.gro \
  -f trajectory.gro \
  -n index.ndx \
  --ref LI \
  --sel Cl \
  --rcut 0.35 \
  --BL \
  -o bond_stats.xvg
```

**Output files:**
- `bond_stats.xvg`: Frame numbers and total bond counts
- `bond_stats_BL.xvg`: Frame numbers and average bond lengths

### Example 4: Center-of-Mass Analysis

Use molecular centers of mass instead of individual atoms:

```bash
python3 gmx_bond.py \
  -s system.gro \
  -f trajectory.gro \
  -n index.ndx \
  --ref CHOL \
  --sel DPPC \
  --rcut 0.8 \
  --ref_mol true \
  --sel_mol true \
  --mass_file masses.txt \
  -o mol_bonds.xvg
```

**Key points:**
- Requires `--mass_file` with atomic masses
- Counts bonds between molecules instead of atoms
- Useful for coarse-grained analysis

### Example 5: Frame Range and Skip

Analyze specific frame range with stride:

```bash
python3 gmx_bond.py \
  -s system.gro \
  -f trajectory.gro \
  -n index.ndx \
  --ref LI \
  --sel Cl \
  --rcut 0.35 \
  --begin 100 \
  --end 1000 \
  --skip 10 \
  -o frame_subset.xvg
```

**Processing:**
- Starts from frame 100
- Stops at frame 1000
- Reads every 10th frame

### Example 6: Debug Mode

Enable detailed debugging output:

```bash
python3 gmx_bond.py \
  -s system.gro \
  -f trajectory.gro \
  -n index.ndx \
  --ref LI \
  --sel Cl \
  --rcut 0.35 \
  --debug \
  -o debug_bonds.xvg
```

**Outputs:**
- Detailed group information in log
- Bond details for first reference atom in first frame
- Index lists and distances for bonded atoms

---

## Mass File Format

Required for center-of-mass calculations (`--ref_mol true` or `--sel_mol true`).

**Format:**
```
[GROUP_NAME]
atom_name = mass_value
atom_name = mass_value
...

[ANOTHER_GROUP]
atom_name = mass_value
```

**Example (masses.txt):**
```
[CHOL]
C = 12.01
H = 1.008
O = 16.00
N = 14.01

[DPPC]
C = 12.01
H = 1.008
O = 16.00
P = 30.97
```

---

## Understanding the Output

### Bond Count Interpretation

**Total bonds** = Sum of all bonds detected in a frame

For a Li-Cl pair with 50 Li atoms and 75 Cl atoms:
- If each Li coordinates 2 Cl on average: 50 × 2 = 100 total Li-Cl bonds
- This is the actual count, not normalized

### Bond Statistics in Log

```
Bond statistics for pair: LI-Cl
  N_ref (LI): 50
  N_sel (Cl): 75
  Mean total bonds per frame: 102.45
  Min total bonds per frame: 98
  Max total bonds per frame: 106
```

**Interpretation:**
- 50 Li atoms are the reference particles
- 75 Cl atoms are the selection particles
- On average, 102.45 Li-Cl bonds per frame
- Ranges from 98 to 106 bonds

### Cutoff Distance Selection

**Tips for choosing rcut:**
- Use radial distribution function (RDF) peaks
- Typical ionic distances: 0.3-0.5 nm
- Hydrogen bonds: 0.15-0.35 nm
- Check: `gmx rdf -f trajectory.gro -s system.tpr -n index.ndx`

---

## Troubleshooting

### Issue: "Index group not found"
**Solution:** Check group names in index file match `--ref` and `--sel` arguments (case-sensitive)

### Issue: "Zero count for reference/selection"
**Causes:**
- Group name misspelled
- Empty group in index file
- Selected atoms outside trajectory bounds

### Issue: "Mass file not parsed correctly"
**Solution:** Verify mass file format matches example above

### Issue: "Reached EOF before frame"
**Solution:** Trajectory has fewer frames than requested. Reduce `--end` value

### Issue: Memory error with large trajectories
**Solution:** Use `--skip` parameter to read fewer frames, or process in smaller chunks

---

## Performance Tips

1. **Use `--skip` parameter** to reduce number of frames processed
2. **Minimize group sizes** if possible (fewer atoms = faster computation)
3. **Use `--begin` and `--end`** to process only relevant frame range
4. **Disable `--BL`** unless bond length statistics are needed
5. **COM calculation** (`--ref_mol true`) is slower due to mass-weighted averaging

**Performance:** ~100-1000 frames per second depending on group sizes and hardware

---

## Citation

If you use this tool in your research, please cite:
- GROMACS: Abraham et al., SoftwareX 1-2 (2015) 19-25
- Minimal Image Convention: Allen & Tildesley, Computer Simulation of Liquids (1987)

---

## Related Scripts

- `gmx_rdf.py` - Radial distribution function analysis
- `gmx_HbondAnalysis.py` - Hydrogen bond detection and analysis
- `gmx_angle.py` - Bond angle analysis

---

## Contact & Support

For issues or questions, refer to the repository issues page or documentation.

---

## Version History

**v1.0 (2026-04-14)**
- Initial release
- Support for atom and COM-based analysis
- Vectorized bond detection
- Total bond count output (not averaged)
- Optional bond length statistics
