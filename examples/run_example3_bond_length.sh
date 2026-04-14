#!/bin/bash
#
# Example 3: Bond length statistics
# Analyze both total bond counts and average bond lengths
#

echo "Running Example 3: Bond Length Statistics"
echo "=========================================="

python3 gmx_bond.py \
  -s examples/ionic_liquid/system.gro \
  -f examples/ionic_liquid/trajectory.gro \
  -n examples/ionic_liquid/index.ndx \
  --ref LI \
  --sel Cl \
  --rcut 0.35 \
  --BL \
  -o output/example3_bond_stats.xvg

echo ""
echo "This generates two output files:"
echo "  1. output/example3_bond_stats.xvg - Total bond counts"
echo "  2. output/example3_bond_stats_BL.xvg - Average bond lengths"
echo ""
echo "Log file: output/example3_bond_stats.log"
echo ""
echo "To visualize both:"
echo "  xmgrace output/example3_bond_stats.xvg output/example3_bond_stats_BL.xvg"
