# `fix pmf/umbrella/z` for LAMMPS

A custom LAMMPS fix for umbrella sampling along the center-of-mass `z` coordinate of a group, for example a polymer approaching a substrate or electrode.

The fix applies a harmonic bias,

```text
U_bias(z) = 0.5*K*(z_COM - z0)^2
F_spring,z = -K*(z_COM - z0)
```

It automatically moves through a sequence of umbrella windows and writes:

- a summary file with mean COM position, mean spring force, and mean bias energy;
- a histogram file that can be analyzed with WHAM/MBAR to obtain a PMF.

## Files

```text
fix_pmf_umbrella_z.cpp
fix_pmf_umbrella_z.h
wham_pmf_umbrella_z.py
```

## Installation

Copy the C++ files into a custom LAMMPS source directory, for example:

```bash
cp fix_pmf_umbrella_z.cpp /path/to/lammps/src/USER-CUSTOM/
cp fix_pmf_umbrella_z.h   /path/to/lammps/src/USER-CUSTOM/
```

Then rebuild LAMMPS.

## Syntax

```lammps
fix ID group-ID pmf/umbrella/z K zstart zstop dz nequil nprod nevery binwidth file prefix
```

Arguments:

```text
K         harmonic spring constant
zstart    first umbrella-window center
zstop     last umbrella-window center
dz        spacing between neighboring windows; use a positive value
nequil    equilibration steps per window
nprod     production steps per window
nevery    sample histogram every this many production steps
binwidth  histogram bin width for z_COM
prefix    output-file prefix
```

The sign of the scan direction is chosen automatically from `zstart` and `zstop`.

## Example

```lammps
compute com polymer com

fix pmf polymer pmf/umbrella/z 2.0 46.0 -46.0 2.0 100000 200000 100 0.25 file pmf
fix_modify pmf energy yes

thermo 1000
thermo_style custom step temp pe etotal c_com[3] f_pmf f_pmf[1] f_pmf[2] f_pmf[3] f_pmf[4] f_pmf[5] f_pmf[6] f_pmf[7] f_pmf[8]
thermo_modify flush yes norm no

run 14100000
```

This scans:

```text
z0 = 46, 44, 42, ..., -44, -46
```

There are 47 windows. Each window runs `100000 + 200000 = 300000` steps, so the total run length is:

```text
47 * 300000 = 14100000 steps
```

## Output files

For `file pmf`, the fix writes:

```text
pmf.summary.dat
pmf.hist.dat
```

### `pmf.summary.dat`

Columns:

```text
window z0 nsamples <z_COM> std_z <F_spring_z> <F_system_z> <U_bias> out_low out_high
```

where

```text
<F_system_z> = -<F_spring_z>
```

This file is useful for restrained mean-force/umbrella-integration analysis.

### `pmf.hist.dat`

Columns:

```text
window z0 K bin_center count
```

This file is used by WHAM/MBAR to reconstruct the unbiased PMF.

## WHAM analysis

For LJ units with `temp = 1`, run:

```bash
python wham_pmf_umbrella_z.py pmf.hist.dat --kbt 1.0 -o pmf_profile.dat
```

For real units at 300 K, run:

```bash
python wham_pmf_umbrella_z.py pmf.hist.dat --units real --temperature 300 -o pmf_profile.dat
```

The output `pmf_profile.dat` contains:

```text
z PMF probability total_counts
```

The PMF is shifted so that its minimum is zero.

## Notes

This fix does not magically make PMF calculations cheap. Umbrella sampling still requires enough overlap between neighboring windows and sufficient equilibration of polymer conformations, ions, charge regulation, and electrode charges.

Useful checks:

```bash
xmgrace -block pmf.summary.dat -bxy 2:4   # <z_COM> vs z0
xmgrace -block pmf.summary.dat -bxy 2:7   # mean system force vs z0
xmgrace pmf_profile.dat                   # WHAM PMF
```

If `<z_COM>` is far from `z0`, increase `K`, reduce `dz`, or increase `nequil`.

If neighboring histograms do not overlap, reduce `dz`, reduce `K`, or sample longer.
