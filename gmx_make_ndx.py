#!/usr/bin/env python3
import argparse


def generate_index(mol_apm_nmol_list, output_file):
    """Generate index file with sequential atoms (old mode - apm is count)."""
    with open(output_file, 'w') as f:
        current_index = 1
        for mol_name, apm, nmol in mol_apm_nmol_list:
            f.write(f"[ {mol_name} ]\n")
            for i in range(nmol):
                for j in range(apm):
                    f.write(f"{current_index}\n")
                    current_index += 1


def generate_index_explicit(mol_name, apm_atoms, nmol, amol, output_file):
    """Generate index file with explicit atoms (new mode - apm is list of atom indices)."""
    with open(output_file, 'w') as f:
        f.write(f"[ {mol_name} ]\n")
        for mol_id in range(nmol):
            for atom in apm_atoms:
                global_atom = mol_id * amol + atom
                f.write(f"{global_atom}\n")


def main():
    parser = argparse.ArgumentParser(description="Generate GROMACS index file for molecules")
    parser.add_argument("-mol", nargs='+', help="List of molecule names")
    parser.add_argument("-apm", nargs='+', type=int, help="Atoms per molecule OR explicit atom indices (depends on -amol)")
    parser.add_argument("-nmol", nargs='+', type=int, help="Number of molecules for each -mol")
    parser.add_argument("-amol", nargs='+', type=int, default=None, help="Total atoms per molecule (activates explicit atoms mode)")
    parser.add_argument("-o", default="index.ndx", help="Output index file (default: index.ndx)")
    args = parser.parse_args()

    # Check if explicit atoms mode is enabled (when -amol is provided)
    if args.amol is not None:
        # Explicit atoms mode: -apm contains explicit atom indices
        if len(args.mol) > 1:
            parser.error("Explicit atoms mode (-amol) currently supports only single -mol")
        if len(args.nmol) != 1:
            parser.error("Explicit atoms mode requires single -nmol value")
        if len(args.amol) != 1:
            parser.error("Explicit atoms mode requires single -amol value")
        
        mol_name = args.mol[0]
        apm_atoms = args.apm  # Now these are explicit atom indices
        nmol = args.nmol[0]
        amol = args.amol[0]
        
        generate_index_explicit(mol_name, apm_atoms, nmol, amol, args.o)
    else:
        # Original mode: -apm contains atoms per molecule counts
        if not (len(args.mol) == len(args.apm) == len(args.nmol)):
            parser.error("-mol, -apm, and -nmol must have the same number of arguments")

        mol_apm_nmol_list = list(zip(args.mol, args.apm, args.nmol))
        generate_index(mol_apm_nmol_list, args.o)


if __name__ == "__main__":
    main()

