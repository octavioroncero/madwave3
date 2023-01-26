dir=../SRC
dir2=$dir/AUXILIARy-CODES
dir3=$dir/bound-lanz
dir4=$dir/photo

pot=../SRC

######################  compiling modules

      mpif77 -c  -O3  $dir/mod_gridYpara_01y2.f
      mpif77 -c  -O3  $dir/mod_pot_01y2.f
      mpif77 -c  -O3  $dir/mod_baseYfunciones_01y2.f
      mpif77 -c  -O3  $dir/mod_Hphi_01y2.f
      mpif77 -c  -O3  $dir/mod_photoini_01y2.f
      mpif77 -c  -O3  $dir/mod_colini_01y2.f
      mpif77 -c  -O3  $dir/mod_absorcion_01y2.f
      mpif77 -c  -O3  $dir/mod_flux_01y2.f
      mpif77 -c  -O3  $dir/mod_coortrans01_02.f
      mpif77 -c  -O3  $dir3/mod_lanczos_01y2.f

######################  mad3.out
rm mad3.out

mpif77 -O3 -o mad3.out  mod_colini_01y2.o  mod_Hphi_01y2.o mod_gridYpara_01y2.o\
       mod_pot_01y2.o mod_baseYfunciones_01y2.o mod_photoini_01y2.o \
       mod_absorcion_01y2.o mod_flux_01y2.o mod_coortrans01_02.o\
       $dir/main_madwave3.f $dir/liboctdyn.f \
       $dir/fit_general.f $dir/dipele_general.f $dir/coupling_general.f\
 -lfftw3
######################  mad3.out
rm bndgrid.out
mpif77 -O3  -o bndgrid.out  mod_Hphi_01y2.o mod_gridYpara_01y2.o\
       mod_pot_01y2.o mod_baseYfunciones_01y2.o mod_photoini_01y2.o\
       mod_absorcion_01y2.o mod_lanczos_01y2.o\
       $dir3/main_boundlanz.f $dir/liboctdyn.f $dir3/liboctdynlanz.f \
       $dir/fit_general.f $dir/dipele_general.f $dir/coupling_general.f \
 -lfftw3

######################  distri.out & distriREAC.out
rm distri.out distriREAC.out
gfortran  -g -C -o distri.out $dir2/distriwvp.f $dir/liboctdyn.f
gfortran  -O3 -o distriREAC.out $dir2/distriREACwvp.f $dir/liboctdyn.f

######################  crp.out & cip.out & cipref.out
rm cip.out crp.out cipref.out
gfortran  -O3 -o crp.out $dir2/CRP-fast.f $dir/liboctdyn.f
gfortran  -O3 -o cip.out $dir2/CIP-fast.f $dir/liboctdyn.f
gfortran  -O3 -o cipref.out $dir2/CIP-OmgRef.f $dir/liboctdyn.f

######################  rate.out & rates2s.out
rm rate.out rates2s.out
gfortran  -O3 -o rate.out $dir2/rateFromSigma.f $dir/liboctdyn.f
gfortran  -O3 -o rates2s.out $dir2/rate-s2s-fromCRP-extrapolation.f
gfortran  -O3 -o inelastic-rates2s.out $dir2/rate-s2s-fromCIP-extrapolation.f

######################  sigma.out 
rm sigma.out
mpif77  -O3 -o sigma.out $dir2/sigmaFromS2prod.f $dir/liboctdyn.f
######################  cheby-spectra.out 
rm cheby-spectra.out
gfortran  -O3 -o cheby-spectra.out $dir4/cheby-spectra.f 

######################  Einstein.out
rm Einstein.out
mpif77 -o Einstein.out -O3  mod_gridYpara_01y2.o mod_baseYfunciones_01y2.o mod_pot_01y2.o mod_Hphi_01y2.o mod_photoini_01y2.o $dir4/Acoefficient-bound-bound.f $dir/liboctdyn.f $dir/fit_general.f $dir/dipele_general.f $dir/coupling_general.f -lfftw3


echo "13 executable codes: mad3.out bndgrid.out Einstein.out distri.out  distriREAC.out  crp.out  cip.out cipref.out inelastic-rates2s.out rate.out rates2s.out sigma.out  cheby-spectra.out"

echo "removing *.o and *.mod"
rm -f *.o *.mod
