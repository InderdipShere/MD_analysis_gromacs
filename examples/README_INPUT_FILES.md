# Sample input files for Qn_order_parameter

This directory shows example input file formats.

## Files Needed

1. **structure.gro** - GROMACS structure file
2. **trajectory.gro** - GROMACS trajectory file
3. **index.ndx** - GROMACS index file
4. **masses.dat** (optional) - Mass file for COM calculation

## How to Generate These Files

### From GROMACS Trajectory (*.xtc)

```bash
# Export structure and trajectory
gmx editconf -f simulation.tpr -o structure.gro
gmx trjconv -f simulation.xtc -s simulation.tpr -o trajectory.gro

# Create or update index file
gmx make_ndx -f structure.gro -o index.ndx
```

### Create Mass File

```bash
# Extract atomic masses from topology
# Option 1: Manual creation (see masses.dat.example below)
# Option 2: Automatic extraction

gmx dump -s simulation.tpr | grep -A 100 "mass" > masses_tmp.txt
```

## File Format Examples

### structure.gro Example

```
Title
  105
    1SOL      O    1   0.126   0.178   0.245
    1SOL      H    2   0.118   0.181   0.264
    1SOL      H    3   0.109   0.174   0.232
    2SOL      O    4   0.325   0.254   0.156
    2SOL      H    5   0.318   0.261   0.142
    2SOL      H    6   0.311   0.258   0.129
    3NA       NA   7   0.500   0.400   0.300
    4CL       CL   8   0.600   0.350   0.280
    ...
   100SOL      O  100   0.815   0.623   0.489
   100SOL      H  101   0.823   0.619   0.472
   100SOL      H  102   0.806   0.628   0.501
   1.0    1.0    1.0
```

### index.ndx Example

```
[ System ]
     1     2     3     4     5     6     7     8     9    10
    11    12    13    14    15    16    17    18    19    20
    ...

[ SOL ]
     1     2     3     4     5     6     7     8     9    10
    11    12    13    14    15    16    17    18    19    20
    21    22    23    24    25    26    27    28    29    30
    ...

[ NA ]
   301

[ CL ]
   302

[ IONS ]
   301   302
```

### masses.dat Example

```
SOL    18.015
O      15.999
H       1.008
NA     22.990
CL     35.453
```

## Requirements

- .gro files must have consistent residue numbering
- Index groups must have valid atom numbers
- Mass file names should match residue names exactly (case-sensitive)
- All files should use consistent units (nm for distances)

## Verification

Check your input files:

```bash
# Check structure
gmx check -f structure.gro

# Check trajectory compatibility
gmx check -f trajectory.gro -s structure.gro

# List available groups in index
grep "\[" index.ndx
```
