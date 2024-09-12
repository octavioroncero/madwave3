MPI_FC=mpif77
FC=mpif77
FFLAGS=-O3 #-DUSE_MPIF -fallow-invalid-boz -fallow-argument-mismatch # uncomment for newer compilers

SRC_DIR=../SRC/
AUX_DIR=$(SRC_DIR)/AUXILIARy-CODES/
LANZ_DIR=$(SRC_DIR)/bound-lanz/
PHOTO_DIR=$(SRC_DIR)/photo/

MODS=$(addprefix $(SRC_DIR), mod_gridYpara_01y2.f mod_pot_01y2.f mod_baseYfunciones_01y2.f mod_Hphi_01y2.f mod_photoini_01y2.f mod_colini_01y2.f mod_absorcion_01y2.f mod_flux_01y2.f mod_coortrans01_02.f) $(LANZ_DIR)/mod_lanczos_01y2.f
MODS_OBJ=$(MODS:%.f=%.o)

main: mad3 bndgrid
aux: distri distriREAC crp cip cipave rate rates2s inelastic-rates2s sigma cheby-spectra einstein

mad3: $(MODS_OBJ) $(SRC_DIR)/main_madwave3.f $(SRC_DIR)/liboctdyn.f $(SRC_DIR)/fit_general.f $(SRC_DIR)/dipele_general.f $(SRC_DIR)/coupling_general.f
	$(MPI_FC) $(FFLAGS) $^ -o mad3.out -lfftw3

bndgrid: $(MODS_OBJ) $(LANZ_DIR)/main_boundlanz.f $(SRC_DIR)/liboctdyn.f $(LANZ_DIR)/liboctdynlanz.f $(SRC_DIR)/fit_general.f $(SRC_DIR)/dipele_general.f $(SRC_DIR)/coupling_general.f
	$(MPI_FC) $(FFLAGS) $^ -o bndgrid.out -lfftw3

distri: $(AUX_DIR)/distriwvp.f $(SRC_DIR)/liboctdyn.f
	$(FC) $(FFLAGS) -C $^ -o distri.out

distriREAC: $(AUX_DIR)/distriREACwvp.f $(SRC_DIR)/liboctdyn.f
	$(FC) $(FFLAGS) $^ -o distriREAC.out

crp: $(AUX_DIR)/CRP-fast.f $(SRC_DIR)/liboctdyn.f
	$(FC) $(FFLAGS) $^ -o crp.out

cip: $(AUX_DIR)/CIP-fast.f $(SRC_DIR)/liboctdyn.f
	$(FC) $(FFLAGS) $^ -o cip.out
	
cipave: $(AUX_DIR)/CIPaverage-fast.f $(SRC_DIR)/liboctdyn.f
	$(FC) $(FFLAGS) $^ -o cipave.out
	
rate: $(AUX_DIR)/rateFromSigma.f  $(SRC_DIR)/liboctdyn.f
	$(FC) $(FFLAGS) $^ -o rate.out

rates2s: $(AUX_DIR)/rate-s2s-fromCRP-extrapolation.f
	$(FC) $(FFLAGS) $^ -o rates2s.out

inelastic-rates2s: $(AUX_DIR)/rate-s2s-fromCIP-extrapolation.f
	$(FC) $(FFLAGS) $^ -o inelastic-rates2s.out

sigma: $(AUX_DIR)/sigmaFromS2prod.f $(SRC_DIR)/liboctdyn.f
	$(FC) $(FFLAGS) $^ -o sigma.out

cheby-spectra: $(PHOTO_DIR)/cheby-spectra.f
	$(FC) $(FFLAGS) $^ -o cheby-spectra.out

einstein: $(PHOTO_DIR)/Acoefficient-bound-bound.f $(MODS_OBJ) $(SRC_DIR)/liboctdyn.f $(SRC_DIR)/fit_general.f $(SRC_DIR)/dipele_general.f $(SRC_DIR)/coupling_general.f
	$(FC) $(FFLAGS) $^ -o Einstein.out -lfftw3


clean:
	rm -f *.mod *.o *.out $(SRC_DIR)*.o $(AUX_DIR)*o $(LANZ_DIR)*.o $(PHOTO_DIR)*.o 

%.o: %.f
	$(MPI_FC) $(FFLAGS) -c $<  -o $@
