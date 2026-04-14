#!/bin/bash
#
# Example 4: Center-of-Mass (COM) based bond analysis
# Analyze bonds between molecular centers of mass
# Useful for lipids, polymers, or coarse-grained systems
#

echo "Running Example 4: Molecular Center-of-Mass Bond Analysis"
echo "=========================================================="

python3 gmx_bond.py \
  -s examples/lipid_membrane/system.gro \
  -f examples/lipid_membrane/trajectory.gro \
  -n examples/lipid_membrane/index.ndx \
  --ref CHOL \
  --sel DPPC \
  --rcut 0.8 \
  --ref_mol true \
  --sel_mol true \
  --mass_file examples/lipid_membrane/masses.txt \
  -o output/example4_molecular_bonds.xvg

echo ""
echo "This analyzes:"
echo "  - Bonds between cholesterol (CHOL) and DPPC molecules"
echo "  - Using molecular centers of mass (COM)"
echo "  - Cutoff distance: 0.8 nm"
echo ""
echo "Output file: output/example4_molecular_bonds.xvg"
echo "Log file: output/example4_molecular_bonds.log"
echo ""
echo "Key features:"
echo "  - Counts bonds between molecule centers (not atoms)"
echo "  - Requires mass file for COM calculation"
echo "  - Useful for coarse-grained or multi-atom group analysis"
