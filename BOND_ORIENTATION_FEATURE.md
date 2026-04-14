# Bond Orientation Analysis Feature

## Overview

`gmx_bond_orientation.py` extends GROMACS trajectory analysis to calculate **bond orientations** with respect to a reference direction within molecules.

### Key Concepts

- **Bond**: Connection between atoms in reference and selection groups (distance ≤ rcut)
- **Reference Direction**: Defined by vector from atom1 to closest atom2 within each molecule (atom2 can be from different molecule)
- **Orientation Angle**: Angle between bond vector (ref→sel) and reference direction vector, in degrees (0-180°)
- **Filtering**: Only bonds from reference atoms in molecules that have atom1 are analyzed

## Command-Line Interface

```bash
python gmx_bond_orientation.py \
  -s <structure.gro> \
  -f <trajectory.gro> \
  -n <index.ndx> \
  -o <output.xvg> \
  --ref <ref_group1> [<ref_group2> ...] \
  --sel <sel_group1> [<sel_group2> ...] \
  --ref_mol <0|1> [<0|1> ...] \
  --sel_mol <0|1> [<0|1> ...] \
  -direction_species <atom1_ref> <atom2_ref> [<atom1_ref2> <atom2_ref2> ...] \
  --direction_mol <0|1> [<0|1> ...] \
  --rcut <rcut1> [<rcut2> ...] \
  [--mass_file <mass_file.txt>] \
  [--begin <frame>] [--end <frame>] [--skip <N>] \
  [-dist <N_bins>] \
  [-debug]
```

## Required Arguments

| Argument | Description |
|----------|-------------|
| `-s, --structure_file` | Initial structure file (.gro) |
| `-f, --traj_file` | Trajectory file (.gro) |
| `-n, --index_file` | Index file (.ndx) with atom groups |
| `--ref` | Name(s) of reference group(s) - atoms where bonds originate |
| `--sel` | Name(s) of selection group(s) - atoms where bonds terminate |
| `-direction_species` | **2 × Npairs atom names** defining direction vectors (see below) |
| `--rcut` | Cutoff distance(s) in nm for bond detection |

## Optional Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `-o, --output_file` | Output file for orientation statistics | `orientation.xvg` |
| `--ref_mol` | Use center-of-mass for ref groups (0/1) | 0 (atom-based) |
| `--sel_mol` | Use center-of-mass for sel groups (0/1) | 0 (atom-based) |
| `--direction_mol` | Use COM for direction atoms (0/1) | 0 (atom-based) |
| `--mass_file` | File with atomic masses (required if using COM) | None |
| `--begin` | First frame to process | 0 |
| `--end` | Last frame to process (-1 = end) | -1 |
| `--skip` | Process every Nth frame | 1 |
| `-dist` | Generate distribution with N bins (0-180°) | None |
| `-debug` | Enable debug output | False |

## Understanding `-direction_species`

The `-direction_species` option specifies atoms that define the reference direction **for each bond pair**.

For each molecule that contains atom1:
- Find the closest atom2 (can be from any molecule, only needs to be within rcut)
- Compute direction vector: atom1 → atom2 (applying MIC)
- Use this direction to calculate orientation angles for all bonds from that molecule's reference atoms

**Format:** Provide **exactly 2 × Npairs** atom names:
- First pair: atoms defining direction for ref[0]-sel[0]
- Second pair: atoms defining direction for ref[1]-sel[1]
- And so on...

**Example:**
```bash
# One bond pair: Ch with ChCl_Cl, direction from ChCl_N to ChCl_C2
--ref Ch --sel ChCl_Cl \
-direction_species ChCl_N ChCl_C2

# Two bond pairs
--ref Ch ChCl_N --sel ChCl_Cl ChCl_Cl \
-direction_species ChCl_N ChCl_C2 ChCl_O ChCl_N
```

## Molecule-based Filtering

**Key behavior:**
- For each molecule, the script extracts atoms matching `-direction_species` atom names
- If both atoms exist in the molecule, a direction vector is computed (atom1 → atom2)
- Bonds are calculated **only for reference atoms from molecules with valid direction vectors**
- Example: If you specify `--ref_mol 1` (COM of molecules), the code calculates for each reference molecule separately

## Output Files

### 1. Main Orientation File (`.xvg`)

**Default:** `orientation.xvg`

Contains:
- **Column 1:** Frame number
- **Columns 2+:** Average orientation angle per frame for each bond pair (in degrees)

**Example output:**
```
# frame Orientation_Ch-ChCl_Cl(deg)
0       45.234567
1       46.123456
2       44.987654
...
```

**Footer comments:** Total numbers of reference and selection particles per pair

### 2. Distance Distribution File (`_dist.xvg`) - Optional

**Generated with:** `-dist N` flag

Contains histogram of orientation angles:
- **Column 1:** Angle (degrees, center of bin)
- **Columns 2+:** Histogram height (normalized) per pair

**Key features:**
- N bins uniformly spaced from 0° to 180°
- Normalized so that integral H(θ)·dθ = **total number of bonds**
- Allows comparison of angle distributions regardless of total bond count
- Log output includes integration check and statistics

### 3. Log File (`.log`)

Comprehensive logging including:
- System information (hostname, platform, Python version)
- Input parameters and configuration
- Group sizes and counts
- Orientation statistics (mean, min, max, std dev per pair)
- Execution timing

## Example Usage

### Basic Usage

```bash
python gmx_bond_orientation.py \
  -s traj.gro \
  -f traj.gro \
  -n index.ndx \
  --ref Ch --sel ChCl_Cl \
  -direction_species ChCl_N ChCl_C2 \
  --rcut 0.35 \
  -o bond_orient.xvg
```

### With COM-based Molecules

```bash
python gmx_bond_orientation.py \
  -s traj.gro \
  -f traj.gro \
  -n index.ndx \
  --ref Ch --sel ChCl_Cl \
  --ref_mol 1 --sel_mol 0 \
  -direction_species ChCl_N ChCl_C2 \
  --direction_mol 1 \
  --mass_file masses.txt \
  --rcut 0.35
```

### With Orientation Distribution

```bash
python gmx_bond_orientation.py \
  -s traj.gro \
  -f traj.gro \
  -n index.ndx \
  --ref Ch --sel ChCl_Cl \
  -direction_species ChCl_N ChCl_C2 \
  --rcut 0.35 \
  -dist 36  # 36 bins = 5° per bin
```

### Multiple Bond Pairs

```bash
python gmx_bond_orientation.py \
  -s traj.gro \
  -f traj.gro \
  -n index.ndx \
  --ref Ch ChCl_N --sel ChCl_Cl ChCl_Cl \
  --ref_mol 1 0 --sel_mol 0 0 \
  -direction_species ChCl_N ChCl_C2 ChCl_O ChCl_N \
  --direction_mol 1 0 \
  --mass_file masses.txt \
  --rcut 0.35 0.40 \
  -dist 36
```

## Mathematical Details

### Orientation Angle Calculation

For each bond detected between reference atom i and selection atom j:

1. **Bond vector:** $\vec{b} = \vec{r}_{sel,j} - \vec{r}_{ref,i}$
2. **Direction vector:** $\vec{d} = \vec{r}_{atom2} - \vec{r}_{atom1}$ (or their COM)
3. **Normalize:** $\hat{b} = \vec{b}/|\vec{b}|$, $\hat{d} = \vec{d}/|\vec{d}|$
4. **Angle:** $\theta = \arccos(\hat{b} \cdot \hat{d})$ (in radians)
5. **Convert to degrees:** $\theta_{deg} = \theta \times 180/\pi$

Result: $\theta \in [0°, 180°]$

### Histogram Normalization

When using `-dist N`:

$$H(\theta) = \frac{\text{count}(\theta \text{ bin})}{\Delta\theta}$$

where:
- $\Delta\theta$ = bin width = $180°/N$ 
- Histogram height is normalized by bin width
- Integral: $\int_0^{180°} H(\theta) \, d\theta = N_{\text{total}}$ (total number of bonds)

## Performance Considerations

- **Vectorized calculations:** Uses NumPy broadcasting for efficient distance and angle computations
- **Per-frame processing:** Trajectory processed frame-by-frame to manage memory
- **Sparse operations:** Only computes angles for detected bonds (within rcut)
- **Typical speed:** ~0.1-1.0 seconds per frame depending on system size

## Troubleshooting

| Problem | Solution |
|---------|----------|
| "Could not find atoms ChCl_N or ChCl_C2" | Check atom names in structure file match `-direction_species` |
| "No angles found" | Check rcut is large enough; verify atoms are bonding |
| Mass file errors | Ensure mass file format matches; use `-debug` to inspect groups |
| Integration check ≠ total bonds | Small numerical variation is normal; should match within ±1% for N_bins ≥ 20 |
| atom2 not found in molecule | atom2 can be from any molecule; check it exists in structure |

## Implementation Notes

- **Minimum image convention (MIC):** Applied for all distance/vector calculations across PBC
- **Direction vector:** Each molecule computes direction from atom1 (in molecule) to closest atom2 (any molecule)
- **Filtering:** Only bonds from reference atoms in molecules that have atom1 are analyzed
- **Frame averaging:** Each frame outputs average angle across all valid bonds
- **Histogram normalization:** Integral equals total number of bonds analyzed
- **Statistics:** Mean, min, max, and standard deviation logged for each pair

## Related Scripts

- [gmx_bond.py](gmx_bond.py) - Basic bond statistics (count, length, distance distribution)
- [gmx_bond_orientation.py](gmx_bond_orientation.py) - This script (orientation analysis)
