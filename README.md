# Madwave3: quantum wave packet program for triatomic systems  

It is a fortran code, parallelized with MPI and openMP,
for the quantum propagation of wave packets describing the dynamics of triatomic systems,
 for treating state-to-state reactive and inelastic scattering::

   $ 01(v,j) + 2 --> 02(v',j') + 1

and photodissociation::

   $ 01-2 +h nu --> 01(v,j) +2
   $            --> 02(v',j') + 1

in several electronic states, which are considered to be
diabatic (diagonal) in the 01+2 reactant channel asymptote, and non-diagonal
for 02+1 product channel, as a general situation.

The user-provided potential must be set for each particular 012 system,
 according to the examples. For photodissociation, the electric-dipole
transition moments between the initial (adiabatic) and the final
electronic states need also to be be provided.

The main program is mad3.out for the wave packet propagation, but
also a bound state program, bndgrid.out, to calculate the initial vibrational
state for photodissociation. There are other auxiliary programs to
calculate the reactive/inelastic cross sections

## Installation

All source programs (in fortran) and examples are downloaded by cloning
the repository::

```
 git clone https://github.com/qmolastro/madwave3
```

that will create the directory madwave3 with 3 sub-directories

```
      BIN  EXAMPLES  SRC
```

For the instalation of the general purpose program::

```
$ cd BIN 
$ source ./colmad3.sh
```
or
```
$ cd BIN 
$ make main aux
```

the main (compiled with `make main`) programs are `mad3.out` and `bndgrid.out`
devoted to the wavepacket propagation and bound state calculations.

A set of auxiliary programs are also generated (with `make aux`) devoted to the analysis and calculations of different properties after the wavepacket propagation. They are independent from the potential used. These programs are:
- distri.out
- distriREAC.out
- crp.out
- cip.out
- cipave.out
- rate.out
- rates2s.out
- inelastic-rates2s.out
- sigma.out
- cheby-spectra.out
- Einstein.out

## Usage
mad3.out (and bndgrid.out) read the potential, fragments wave functions and electric dipole moments in

```
../pot  ../func ../dip ../bnd  (the  last two  in the case of photodissociation)
```

in which the user-provided potential program writes the required information.

## Examples
An example for the H+HD collision using an in-house PES can be found in the directory
[EXAMPLES/H+DH-v0j0](EXAMPLES/HD%2BH_collision). Please follow the instructions in the README file inside this folder for a detailed description of the input for each calculation and the steps to reproduce them.
