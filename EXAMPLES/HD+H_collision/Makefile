MPI_FC=mpif77
FC=mpif77
FFLAGS=-O3 #-DUSE_MPIF -fallow-invalid-boz -fallow-argument-mismatch # uncomment for newer compilers

SRC_DIR=../../SRC/

FLIB=-lfftw3

MODS=$(addprefix $(SRC_DIR), mod_gridYpara_01y2.f mod_pot_01y2.f mod_baseYfunciones_01y2.f mod_Hphi_01y2.f mod_photoini_01y2.f mod_colini_01y2.f liboctdyn.f dipele_general.f coupling_general.f)

POT_DIR=pes/
POT=$(addprefix $(POT_DIR), FFmodules.f90  FFpots.f90  FFroutines.f  pothead.f)

pot: $(MODS) $(POT) $(SRC_DIR)main_potini.f
	$(MPI_FC) $(FFLAGS) $^ -o pot.out $(FLIB)

clean:
	rm -f pot.out $(SRC_DIR)/*.o $(SRC_DIR)/*.mod
