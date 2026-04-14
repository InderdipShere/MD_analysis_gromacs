#!/bin/bash
#
# Example 2: Multiple bond types
# Analyze several different bond types simultaneously
#

echo "Running Example 2: Multiple Bond Types Analysis"
echo "================================================"

python3 gmx_bond.py \
  -s examples/mixed_system/system.gro \
  -f examples/mixed_system/trajectory.gro \
  -n examples/mixed_system/index.ndx \
  --ref LI Cl Na \
  --sel Cl Na LI \
  --rcut 0.35 0.40 0.38 \
  -o output/example2_multi_bonds.xvg

echo ""
echo "This analyzes three bond types:"
echo "  1. LI-Cl with cutoff 0.35 nm"
echo "  2. Cl-Na with cutoff 0.40 nm"
echo "  3. Na-LI with cutoff 0.38 nm"
echo ""
echo "Output file: output/example2_multi_bonds.xvg"
echo "Log file: output/example2_multi_bonds.log"
