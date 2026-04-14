#!/bin/bash
#
# Example 6: Debug mode
# Enable detailed output for troubleshooting and validation
#

echo "Running Example 6: Debug Mode Analysis"
echo "======================================"

python3 gmx_bond.py \
  -s examples/ionic_liquid/system.gro \
  -f examples/ionic_liquid/trajectory.gro \
  -n examples/ionic_liquid/index.ndx \
  --ref LI \
  --sel Cl \
  --rcut 0.35 \
  --debug \
  -o output/example6_debug_bonds.xvg

echo ""
echo "Debug output includes:"
echo "  - Detailed group atom indices"
echo "  - Molecule composition (if using COM)"
echo "  - Individual bond information for first reference atom"
echo "  - Bond distances and bonded atom/molecule indices"
echo ""
echo "Output file: output/example6_debug_bonds.xvg"
echo "Log file: output/example6_debug_bonds.log (contains detailed debug info)"
echo ""
echo "Use this to:"
echo "  - Verify correct atom group selection"
echo "  - Check bond detection for suspicious results"
echo "  - Validate cutoff distance choices"
echo "  - Troubleshoot unusual output patterns"
