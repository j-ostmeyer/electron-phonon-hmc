# HMC simulations of electron-phonon systems

by Johann Ostmeyer,\
contributions by Pavel Buividovich

An implementation of Hybrid Monte Carlo (HMC) simulations for organic molecular semiconductors and similar electron-phonon systems. The code uses highly efficient Fourier acceleration (FA). Transient localisation (TL) simulations are supported, too.

This repository contains the code required to reproduce the results presented in *"First-principle quantum Monte-Carlo study of charge carrier mobility in organic molecular semiconductors"*, [arXiv:2312.14914 [cond-mat.mtrl-sci]](https://arxiv.org/abs/2312.14914) by Johann Ostmeyer, Tahereh Nematiaram, Alessandro Troisi, Pavel Buividovich.

For questions concerning the code contact [ostmeyer@hiskp.uni-bonn.de](mailto:ostmeyer@hiskp.uni-bonn.de).

## Implementation

The entire simulations are implemented in `C` (located in the `simulation` directory).

### Compile
The code depends on `Lapacke` and `FFTW3`. Make sure you have recent installations available before you proceed.

After adjusting the Makefile according to your needs (possibly switching to another compiler), simply type `make` to compile all the `C` files.

### Run Code
Adjust the `input.txt` file (current default simulates Rubrene at room temperature using HMC with FA) according to the next subsection. Alternatively, create a new file with the relevant input parameters. Then, for HMC simulations, execute:
```
./organic_main input.txt results.csv corr.csv
```

Scalar results will be written to `results.csv` (or any other file provided as a second argument). The Euclidean correlator is written to `corr.csv`. See `organic_hmc.c` and comments therein for details.

For TL simulations execute:
```
./organic_main input.txt results.csv greens.bin
./organic_conv greens.bin smear.txt greens_smeared.csv corr_smeared.csv
```

Again, scalar results are written to `results.csv`. In addition, binned Greens function values will be written to `greens.bin` in binary format. They are analysed by the `organic_conv` execution where, following the instructions in `smear.txt`, the Greens function is convoluted with a smearing kernel for the optical conductivity in `greens_smeared.csv` and transformed to Euclidean time in `corr_smeared.csv`. See `organic_conv.c` and comments therein for details.

### Input file
The input file is structured as follows (comments after the `#` are for this README only and should not be added into the actual input file):
```
L1 L2		# spatial geometry in 2D
Nt		# number of (Euclidean) time slices for HMC / number of Greens function bins for TL
beta		# inverse temperature
J1 J2 J3	# J_a hopping amplitudes on triangular lattice
l1 l2 l3	# lambda_a hopping amplitude fluctuations
mu0		# reduced chemical potential
omega0		# phonon frequency
Nmd traj	# molecular dynamics steps, trajectory length in HMC (traj=0 sets default, optimal for FA) / ignored in TL
therm		# thermalisation steps / ignored in TL
N freq		# number of simulated trajectories, measurement frequency
a1 a2		# 1st unit cell vector in 2D
b1 b2		# 2nd unit cell vector in 2D
flavors		# fermionic flavors
# more flags... (e.g. STATIC_DISORDER for TL simulations)
# see organic_flags.c for options
```

For TL simulations an additional `smear.txt` file is needed:
```
beta Nt		# inverse temperature, number of (Euclidean) time slices for correlator output
w1		# any number of kernel widths for the smearing
w2
...
```

## Data

The data used in the paper will gladly be provided upon request.

## Stable Releases

`v1.0.1` initial publication after the preprint, compatible with arXiv version v2.
