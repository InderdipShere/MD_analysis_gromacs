# Bond Distance Distribution Feature

## Overview

The `-dist N` flag has been successfully added to `gmx_bond.py` to generate **bond distance distribution histograms**.

## New Feature

### Command-Line Argument

```bash
-dist N, --dist N
```

Where `N` is the number of bins for the histogram.

### Output File

- **File name:** `{output_prefix}_dist.xvg`
- **Format:** Multi-column XVG format
- **X-axis:** Distance (nm) from 0 to max(rcut)
- **Y-axis:** Histogram counts (frequency)
- **Columns:** One column per bond pair plus distance bin centers

### What It Does

1. **Tracks all bonded pair distances** during trajectory processing
2. **Creates binned histograms** with N equal-width bins from 0 to max(rcut)
3. **Outputs histogram data** to `{output_file}_dist.xvg`
4. **Provides statistics** in the log file (min, max, mean, median distances)

---

## Usage Examples

### Example 1: Basic Usage

```bash
python3 gmx_bond.py \
  -s system.gro \
  -f trajectory.gro \
  -n index.ndx \
  --ref LI \
  --sel Cl \
  --rcut 0.35 \
  --dist 50 \
  -o bonds.xvg
```

**Output files:**
- `bonds.xvg` - Total bond counts per frame
- `bonds_dist.xvg` - Distance distribution histogram (50 bins)

### Example 2: Combined with Bond Length Statistics

```bash
python3 gmx_bond.py \
  -s system.gro \
  -f trajectory.gro \
  -n index.ndx \
  --ref LI \
  --sel Cl \
  --rcut 0.35 \
  --BL \
  --dist 50 \
  -o bonds.xvg
```

**Output files:**
- `bonds.xvg` - Total bond counts per frame
- `bonds_BL.xvg` - Average bond lengths per frame
- `bonds_dist.xvg` - Distance distribution histogram

### Example 3: Multiple Bond Types with Distribution

```bash
python3 gmx_bond.py \
  -s system.gro \
  -f trajectory.gro \
  -n index.ndx \
  --ref LI Cl Na \
  --sel Cl Na LI \
  --rcut 0.35 0.40 0.38 \
  --dist 100 \
  -o multi_bonds.xvg
```

**Output file:**
- `multi_bonds.xvg` - Bond counts for all pairs
- `multi_bonds_dist.xvg` - Distance distributions with 3 columns (one per pair)

---

## Output File Format

### bonds_dist.xvg

```
# distance(nm) hist_LI-Cl hist_Cl-Na hist_Na-LI
0.000000   450  200  100
0.007500   520  250  120
0.015000   480  280  110
...
0.342500    10    5    2
# Total number of Ref:
# LI: 50, Cl: 75, Na: 60
# Total number of sel:
# Cl: 75, Na: 60, LI: 50
```

**Interpretation:**
- **Column 1:** Distance bin centers (nm)
- **Column 2+:** Histogram counts for each pair
- **Comments:** Total number of reference and selection atoms/molecules

---

## Log File Output

When using `-dist`, the log file includes:

```
--- Bond Distance Distribution ---
Distance distribution for pair LI-Cl:
  Total bonded pairs: 5234
  Min distance: 0.251000 nm
  Max distance: 0.349500 nm
  Mean distance: 0.302345 nm
  Median distance: 0.301200 nm
  Histogram bins: 50, range: 0.0-0.35 nm
Distance distribution data written to bonds_dist.xvg
```

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `-dist N` | None | Number of bins for histogram (optional) |
| X-axis range | 0 to max(rcut) | Automatically determined |
| Y-axis | Counts | Histogram frequency |

---

## How It Works

### Distance Tracking

During bond detection:
1. All distances between reference and selection positions are calculated
2. Distances of bonded pairs (≤ rcut) are stored
3. Distances accumulate across all frames

### Histogram Generation

1. X-axis range: 0.0 to max(rcut) nm
2. Number of equal-width bins: N (user specified)
3. Bin width: max(rcut) / N
4. Histogram counts the frequency in each bin

### Multi-pair Support

- Each pair gets its own histogram column
- Different pairs can have different rcut values
- All pairs use the same bin edges (0 to max(rcut))
- Distances beyond the global max(rcut) are not included

---

## Analysis Tips

### Interpreting the Histogram

1. **Peak location:** Shows most common bond distance
2. **Distribution shape:** 
   - Sharp peak → rigid bonds with tight coordination
   - Broad distribution → flexible bonds, varying coordination geometries
3. **Multiple peaks:** Possible multiple coordination shell structures

### Visualization

Use with Grace, Gnuplot, or Python:

```bash
# Grace/xmgrace
xmgrace bonds_dist.xvg

# Python matplotlib
python3 << 'EOF'
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('bonds_dist.xvg', comments='#')
distances = data[:, 0]

# Plot histograms
for i in range(1, data.shape[1]):
    plt.plot(distances, data[:, i], label=f'Pair {i}')

plt.xlabel('Distance (nm)')
plt.ylabel('Count')
plt.legend()
plt.savefig('bond_distribution.png')
EOF
```

---

## Combined Analysis Workflow

```bash
# Step 1: Run full analysis with all features
python3 gmx_bond.py \
  -s system.gro \
  -f trajectory.gro \
  -n index.ndx \
  --ref LI \
  --sel Cl \
  --rcut 0.35 \
  --BL           # Get average bond lengths
  --dist 50      # Get distance distribution
  -o analysis.xvg

# Step 2: Examine outputs
# analysis.xvg          - Total bond counts per frame
# analysis_BL.xvg       - Average bond lengths per frame
# analysis_dist.xvg     - Distance distribution histogram

# Step 3: Visualize
xmgrace analysis.xvg analysis_BL.xvg analysis_dist.xvg
```

---

## Performance Considerations

- **Memory:** Stores all bonded pair distances (can be large for long trajectories)
- **Time:** Minimal overhead (no extra distance calculations)
- **Disk:** Additional output file of size ~N × n_pairs bytes

**Optimization:**
- Use `--skip` parameter to reduce memory usage on large trajectories
- Use `--begin` and `--end` to analyze specific regions

---

## Common Use Cases

### 1. Validate Cutoff Distance

Check if your chosen rcut captures the first coordination shell:

```bash
python3 gmx_bond.py ... --dist 100
# Look at distribution peak in bonds_dist.xvg
# Should see clear separation from 0
```

### 2. Compare Bond Distributions Across Simulations

```bash
# Run for multiple temperatures/conditions
for TEMP in 300 350 400; do
  python3 gmx_bond.py ... -dist 50 -o T${TEMP}_bonds.xvg
done

# Compare: T300_bonds_dist.xvg vs T350_bonds_dist.xvg vs T400_bonds_dist.xvg
```

### 3. Identify Coordination Shell Transitions

Sharp changes in the distribution shape can indicate:
- Phase transitions
- Structural reorganization
- Coordination environment changes

---

## Technical Notes

- Bin edges range from exactly 0.0 to max(rcut)
- Histogram counts include bonds at exact bin edges
- Empty bins appear as 0 in output
- Distances are measured with periodic boundary conditions (MIC)

---

## Integration with Existing Features

The `-dist` flag works seamlessly with:
- ✅ Atom-based analysis (--ref_mol false)
- ✅ Center-of-mass analysis (--ref_mol true)
- ✅ Multiple bond types (multiple --ref/--sel pairs)
- ✅ Different cutoffs per pair (--rcut val1 val2 ...)
- ✅ Bond length statistics (--BL flag)
- ✅ Frame range limiting (--begin/--end)
- ✅ Frame skipping (--skip)
- ✅ Debug mode (--debug)

---

## Troubleshooting

### Issue: Empty histogram (all zeros)

**Cause:** No bonds detected
- Check cutoff distance is appropriate: `gmx rdf`
- Verify group selection: use `--debug`
- Check frame range: are bonds present? (see bonds.xvg)

### Issue: Histogram looks strange

**Cause:** Bin width may be inappropriate
- Use more bins for lower resolution
- Use fewer bins for better statistics
- Increase frame statistics with `--skip 1` (process all frames)

### Issue: Memory error

**Cause:** Too many bonded pairs being tracked
- Use `--skip N` to reduce number of frames
- Use `--begin/--end` to limit trajectory range
- Use `-dist` with fewer bins (if reasonable)

---

## Files Modified

- `gmx_bond.py`: Added `-dist` argument, bond distance tracking, histogram generation function

## Backward Compatibility

✅ Fully backward compatible
- Existing scripts work without changes
- `-dist` flag is purely optional
- No changes to existing output files

---

## Example Output Summary

```
BOND ANALYSIS FROM TRAJECTORY
...
--- Calculation ---
Computing bond statistics using vectorized distance calculations...
Total frames processed: 1000

--- Output ---
Bond count data written to bonds.xvg

--- Bond Distance Distribution ---
Distance distribution for pair LI-Cl:
  Total bonded pairs: 5234
  Min distance: 0.251000 nm
  Max distance: 0.349500 nm
  Mean distance: 0.302345 nm
  Median distance: 0.301200 nm
  Histogram bins: 50, range: 0.0-0.35 nm
Distance distribution data written to bonds_dist.xvg

EXECUTION SUMMARY
Total execution time: 45.23 seconds
Time per frame: 0.0452 seconds
```

---

## Next Steps

1. Use `xmgrace` or `matplotlib` to visualize the histogram
2. Compare distributions across different simulations
3. Use statistical analysis to extract quantitative insights
4. Identify anomalous bond formation patterns

Enjoy your bond analysis! 🎉
