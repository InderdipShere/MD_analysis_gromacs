#!/bin/bash
#
# Master Example Runner
# Demonstrates all available features of gmx_bond.py
#
# Usage: ./run_all_examples.sh
#

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUTPUT_DIR="$SCRIPT_DIR/output"

echo ""
echo "╔════════════════════════════════════════════════════════════════╗"
echo "║         GROMACS Bond Analysis Tool - Example Suite             ║"
echo "║                  All Examples Demonstration                    ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

echo "The following examples demonstrate different use cases:"
echo ""
echo "1. Basic Analysis       - Single bond type (Atom-based)"
echo "2. Multiple Bonds       - Several bond types simultaneously"
echo "3. Bond Length Stats    - Include average bond length output"
echo "4. COM Analysis         - Molecular center-of-mass bonds"
echo "5. Frame Range/Skip     - Process trajectory subsets"
echo "6. Debug Mode           - Detailed troubleshooting output"
echo ""

read -p "Which example would you like to run? (1-6, or 'all'): " choice

case $choice in
  1)
    echo ""
    echo "Running Example 1: Basic Li-Cl Bond Analysis..."
    bash "$SCRIPT_DIR/examples/run_example1_basic.sh"
    ;;
  2)
    echo ""
    echo "Running Example 2: Multiple Bond Types..."
    bash "$SCRIPT_DIR/examples/run_example2_multi_bonds.sh"
    ;;
  3)
    echo ""
    echo "Running Example 3: Bond Length Statistics..."
    bash "$SCRIPT_DIR/examples/run_example3_bond_length.sh"
    ;;
  4)
    echo ""
    echo "Running Example 4: COM-Based Analysis..."
    bash "$SCRIPT_DIR/examples/run_example4_com_analysis.sh"
    ;;
  5)
    echo ""
    echo "Running Example 5: Frame Range & Skip..."
    bash "$SCRIPT_DIR/examples/run_example5_frame_range.sh"
    ;;
  6)
    echo ""
    echo "Running Example 6: Debug Mode..."
    bash "$SCRIPT_DIR/examples/run_example6_debug.sh"
    ;;
  all)
    echo ""
    echo "Running all examples..."
    for i in {1..6}; do
      echo ""
      echo "═════════════════════════════════════════"
      echo "Example $i"
      echo "═════════════════════════════════════════"
      bash "$SCRIPT_DIR/examples/run_example${i}_*.sh" 2>&1 | head -20
      echo ""
    done
    ;;
  *)
    echo "Invalid choice. Please select 1-6 or 'all'."
    exit 1
    ;;
esac

echo ""
echo "═════════════════════════════════════════"
echo "Analysis complete!"
echo "Output files are in: $OUTPUT_DIR"
echo "═════════════════════════════════════════"
echo ""
