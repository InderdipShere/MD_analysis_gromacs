#!/bin/bash
#
# Example 1: Basic Li-Cl bond analysis
# Simple single bond type analysis for an ionic liquid
#

echo "Running Example 1: Basic Li-Cl Bond Analysis"
echo "=============================================="

python3 gmx_bond.py \
  -s examples/ionic_liquid/system.gro \
  -f examples/ionic_liquid/trajectory.gro \
  -n examples/ionic_liquid/index.ndx \
  --ref LI \
  --sel Cl \
  --rcut 0.35 \
  -o output/example1_li_cl_bonds.xvg

echo ""
echo "Output file: output/example1_li_cl_bonds.xvg"
echo "Log file: output/example1_li_cl_bonds.log"
echo ""
echo "To visualize the results:"
echo "  xmgrace output/example1_li_cl_bonds.xvg"
