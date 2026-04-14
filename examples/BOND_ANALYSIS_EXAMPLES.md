# GROMACS Bond Analysis Examples

This directory contains practical examples of how to use `gmx_bond.py` for different types of bond analysis.

## Quick Start

### Option 1: Run Interactive Menu
```bash
bash run_all_examples.sh
```
This opens an interactive menu where you can select which example to run.

### Option 2: Run Individual Example
```bash
bash run_example1_basic.sh
bash run_example2_multi_bonds.sh
bash run_example3_bond_length.sh
bash run_example4_com_analysis.sh
bash run_example5_frame_range.sh
bash run_example6_debug.sh
```

## Examples Overview

### Example 1: Basic Atom-Based Bond Analysis
**File:** `run_example1_basic.sh`

**What it does:**
- Analyzes Li-Cl bonds in an ionic liquid
- Uses individual atoms as reference points (not COM)
- Single bond type analysis
- Outputs total bond counts per frame

**Command:**
```bash
python3 gmx_bond.py \
  -s system.gro -f trajectory.gro -n index.ndx \
  --ref LI --sel Cl \
  --rcut 0.35 \
  -o li_cl_bonds.xvg
```

**Best for:**
- Simple ion pairing analysis
- Small systems
- Initial exploration of bond formation

---

### Example 2: Multiple Bond Types
**File:** `run_example2_multi_bonds.sh`

**What it does:**
- Analyzes three different bond types simultaneously
- Different cutoff distances for each pair
- Shows how to process multiple interactions at once

**Command:**
```bash
python3 gmx_bond.py \
  -s system.gro -f trajectory.gro -n index.ndx \
  --ref LI Cl Na \
  --sel Cl Na LI \
  --rcut 0.35 0.40 0.38 \
  -o multi_bonds.xvg
```

**Pairs analyzed:**
- LI-Cl (0.35 nm)
- Cl-Na (0.40 nm)
- Na-LI (0.38 nm)

**Best for:**
- Complex ionic systems with multiple species
- Comparing different bond types
- Screening multiple interactions

---

### Example 3: Bond Length Statistics
**File:** `run_example3_bond_length.sh`

**What it does:**
- Calculates both bond counts AND average bond lengths
- Produces two output files:
  - `bonds.xvg`: Total bond counts
  - `bonds_BL.xvg`: Average bond lengths

**Command:**
```bash
python3 gmx_bond.py \
  -s system.gro -f trajectory.gro -n index.ndx \
  --ref LI --sel Cl \
  --rcut 0.35 \
  --BL \
  -o bond_stats.xvg
```

**Output:**
- Frame numbers
- Total number of bonds per frame
- Average bond length per frame

**Best for:**
- Structural characterization
- Investigating bond strength/coordination
- Combined coordination and distance analysis

---

### Example 4: Molecular Center-of-Mass Analysis
**File:** `run_example4_com_analysis.sh`

**What it does:**
- Analyzes bonds between molecular centers of mass (not individual atoms)
- Useful for lipid membranes, polymers, coarse-grained systems
- Requires atomic masses for COM calculation

**Command:**
```bash
python3 gmx_bond.py \
  -s system.gro -f trajectory.gro -n index.ndx \
  --ref CHOL --sel DPPC \
  --rcut 0.8 \
  --ref_mol true \
  --sel_mol true \
  --mass_file masses.txt \
  -o molecular_bonds.xvg
```

**Key features:**
- `--ref_mol true`: Use COM for reference molecules
- `--sel_mol true`: Use COM for selection molecules
- `--mass_file`: Atomic masses for COM calculation

**Best for:**
- Lipid-lipid interactions
- Protein-ligand contacts
- Polymer bridging
- Coarse-grained analysis

---

### Example 5: Frame Range and Skip
**File:** `run_example5_frame_range.sh`

**What it does:**
- Processes only a subset of the trajectory
- Skips frames for faster analysis or memory efficiency
- Useful for large production trajectories

**Command:**
```bash
python3 gmx_bond.py \
  -s system.gro -f trajectory.gro -n index.ndx \
  --ref LI --sel Cl \
  --rcut 0.35 \
  --begin 100 \
  --end 1000 \
  --skip 10 \
  -o frame_subset.xvg
```

**Parameters:**
- `--begin 100`: Start at frame 100
- `--end 1000`: Stop at frame 1000
- `--skip 10`: Read every 10th frame (→ ~90 frames processed)

**Best for:**
- Quick analysis of equilibration region
- Processing very long trajectories
- Memory-constrained systems
- Testing analysis parameters

---

### Example 6: Debug Mode
**File:** `run_example6_debug.sh`

**What it does:**
- Enables verbose output for troubleshooting
- Shows detailed atom/molecule lists
- Individual bond information for validation

**Command:**
```bash
python3 gmx_bond.py \
  -s system.gro -f trajectory.gro -n index.ndx \
  --ref LI --sel Cl \
  --rcut 0.35 \
  --debug \
  -o debug_bonds.xvg
```

**Debug output includes:**
- Full list of atom indices in each group
- Molecule information (if using COM)
- Bond details for first reference atom
- Distance values for validation

**Best for:**
- Verifying correct group selection
- Checking for suspicious bonding patterns
- Validating cutoff distances
- Troubleshooting analysis issues

---

## Input Files

### Index File (.ndx)
Create with GROMACS:
```bash
gmx make_ndx -f system.gro -o index.ndx
```

Example structure:
```
[LI]
 1    2    3    4    5

[Cl]
10   11   12   13   14
```

### Mass File (.txt)
Optional, required for COM calculations. See `mass_file_template.txt` for format.

### Structure File (.gro)
GROMACS structure file with atom information.

### Trajectory File (.gro)
GROMACS trajectory in .gro format.

---

## Output Interpretation

### Bond Count Output
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

**Interpretation:**
- `Total_bonds_LI-Cl = 100` means 100 Li-Cl pairs detected in that frame
- **Not normalized** - actual pair count
- Comment lines show reference and selection group sizes

### Bond Length Output
```
# frame BL_LI-Cl(nm)
0      0.280000
1      0.285000
2      0.278000
```

**Note:** Only appears if `--BL` flag is used

---

## Common Workflow

1. **Prepare index file:**
   ```bash
   gmx make_ndx -f system.gro -o index.ndx
   ```

2. **Convert trajectory if needed:**
   ```bash
   gmx trjconv -f traj.xtc -s system.tpr -o traj.gro
   ```

3. **Determine appropriate cutoff:**
   ```bash
   gmx rdf -f traj.gro -s system.tpr -n index.ndx
   ```
   (Check RDF peaks to set `--rcut`)

4. **Run bond analysis:**
   ```bash
   python3 gmx_bond.py -s system.gro -f traj.gro -n index.ndx \
     --ref LI --sel Cl --rcut 0.35 -o bonds.xvg
   ```

5. **Visualize results:**
   ```bash
   xmgrace bonds.xvg
   # or use Grace/plotting software of choice
   ```

---

## Tips & Tricks

### Quick Test
Use `--begin`, `--end`, and `--skip` for quick testing:
```bash
python3 gmx_bond.py [...] --begin 0 --end 50 --skip 5
# Only processes frames 0-50, reading every 5th frame
```

### Determine Cutoff Distance
Always check RDF first:
```bash
gmx rdf -f trajectory.gro -s system.tpr -n index.ndx
```
RDF peaks typically indicate optimal cutoff distances.

### Memory Efficiency
For large trajectories, use stride:
```bash
# Instead of processing all frames:
python3 gmx_bond.py [...] --skip 10  # Every 10th frame
```

### Validate Results
Use debug mode to check first frame:
```bash
python3 gmx_bond.py [...] --debug --begin 0 --end 1
```
This processes only frame 0 with detailed output.

---

## Troubleshooting

### Index group not found
- Check spelling (case-sensitive)
- Verify group exists: `gmx make_ndx -f system.gro -n index.ndx`

### Empty output
- Verify trajectories have frames
- Check cutoff distance isn't too small
- Use `--debug` to see bond details

### Memory errors
- Use `--skip` to reduce frames
- Use `--begin` and `--end` to limit range
- Process in smaller chunks

### Unexpected zero counts
- Verify groups are in trajectory
- Check atom numbers match
- Use `--debug` mode to validate

---

## Next Steps

- Read `BOND_ANALYSIS_GUIDE.md` for detailed documentation
- Experiment with different cutoff distances
- Compare multiple bond types
- Integrate with your analysis workflow

---

## Related Tools

- `gmx_rdf.py` - Radial distribution analysis
- `gmx_HbondAnalysis.py` - Hydrogen bond tracking
- `gmx_angle.py` - Bond angle analysis

Good luck with your analysis!
