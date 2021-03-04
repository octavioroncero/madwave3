====
Madwave3: quantum wavepacket of triatomic systems  
====

It is a fortran77 code, parallelized with MPI and openMP,
for the quantum propagation wave packets describing dynamics in triatomic systems,
 for treating state-to-state reactive and inelastic scattering::

   $ 01(v,j) + 2 --> 02(v',j') + 1

and photodissociation::

   $ 01-2 +h nu --> 01(v,j) +2
   $            --> 02(v',j') + 1

in several electronic states, which are considered to be
diabatic (diagonal) in the 01+2 reactant channel, and non-diagonal
for 02+1 product channel as a general situation.

The user-provided potential must be set for each particular 012 system,
 according to the examples. For photodissociation, the electric-dipole
transition moments between the initial (adiabatic) and the final
electronic states need also to be be provided.

The main program is mad3.out for the wave packet propagation, but
also a bound state program, bndgrid.out, to calculate the initial vibrational
state for photodissociation. There are other auxiliary programs to
calculate the reactive/inelastic cross sections

Installation
============

All source programs (in fortran) and examples are downloaded by cloning
the repository::

$ git clone https://github.com/octavioroncero/madwave3

that will create the directory madwave3 with 4 sub-directories

      BIN  EXAMPLES  PES  SRC

For the instalation of the general purpose program::

$ cd BIN
$ source ./colmad3.sh

and it will create 10 executables

bndgrid.out  cheby-spectra.out  cip.out  colmad3.sh  crp.out  distri.out  distriREAC.out  mad3.out  rate.out  rates2s.out  sigma.out  

which are independent on the potential used. 

mad3.out (and bndgrid.out) read the potential, fragments wave functions and electric dipole moments in

../pot  ../func ../dip ../bnd  (the 2 last in the case of photodissociation)

in which the user-provided potential program write the required information.

An example (for H+HD  using the BKMP2  PES) can be found in directory
EXAMPLES/H+DH-v0j0::

$ ./colpot.sh
$ ./pot.out

will generate that information.
colpot.sh is a shell  that compile the potential to generate pot.out,
providing an example of adapting a external potential

It uses the data in "input.dat" organized in different namelist
which is also used by mad3.out code to calculate state-to-state reation
probabilities for each partial wave (total angular momenbtum J)


