#!/bin/bash
#
# Example 5: Frame range and stride analysis
# Process only a subset of the trajectory with frame skipping
# Useful for very large trajectories or quick analysis
#

echo "Running Example 5: Frame Range & Skip Analysis"
echo "=============================================="

python3 gmx_bond.py \
  -s examples/ionic_liquid/system.gro \
  -f examples/ionic_liquid/trajectory.gro \
  -n examples/ionic_liquid/index.ndx \
  --ref LI \
  --sel Cl \
  --rcut 0.35 \
  --begin 100 \
  --end 1000 \
  --skip 10 \
  -o output/example5_frame_subset.xvg

echo ""
echo "Processing parameters:"
echo "  Start frame: 100"
echo "  End frame: 1000"
echo "  Frame skip: 10 (reads every 10th frame)"
echo "  Total frames processed: ~90"
echo ""
echo "Output file: output/example5_frame_subset.xvg"
echo "Log file: output/example5_frame_subset.log"
echo ""
echo "Use this for:"
echo "  - Quick testing on large trajectories"
echo "  - Processing specific equilibration/production regions"
echo "  - Memory-efficient analysis of long simulations"
