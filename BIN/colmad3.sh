dir=/work/octavio/madwave3/photoBC/SRC
dir2=$dir/AUXILIARy-CODES
dir3=$dir/bound-lanz
dir4=$dir/photo

pot=../SRC

######################  compiling modules

      mpif77 -c  -g -C -fcheck=all  $dir/mod_gridYpara_01y2.f
      mpif77 -c  -g -C -fcheck=all  $dir/mod_pot_01y2.f
      mpif77 -c  -g -C -fcheck=all  $dir/mod_baseYfunciones_01y2.f
      mpif77 -c  -g -C -fcheck=all  $dir/mod_Hphi_01y2.f
      mpif77 -c  -g -C -fcheck=all  $dir/mod_photoini_01y2.f
      mpif77 -c  -g -C -fcheck=all  $dir/mod_colini_01y2.f
      mpif77 -c  -g -C -fcheck=all  $dir/mod_absorcion_01y2.f
      mpif77 -c  -g -C -fcheck=all  $dir/mod_flux_01y2.f
      mpif77 -c  -g -C -fcheck=all  $dir/mod_coortrans01_02.f
      mpif77 -c  -g -C -fcheck=all  $dir3/mod_lanczos_01y2.f
      mpif77 -c  -g -C -fcheck=all  $dir/mod_PsiE_01y2.f
######################  mad3.out
rm mad3.out

mpif77 -g -C -fcheck=all -o mad3.out  mod_colini_01y2.o  mod_Hphi_01y2.o mod_gridYpara_01y2.o\
       mod_pot_01y2.o mod_baseYfunciones_01y2.o mod_photoini_01y2.o \
       mod_absorcion_01y2.o mod_flux_01y2.o mod_coortrans01_02.o mod_PsiE_01y2.o\
       $dir/main_madwave3.f $dir/liboctdyn.f \
       $dir/fit_general.f $dir/dipele_general.f $dir/coupling_general.f\
 -lfftw3
######################  mad3.out
rm bndgrid.out 
mpif77 -g -C -fcheck=all  -o bndgrid.out  mod_Hphi_01y2.o mod_gridYpara_01y2.o\
       mod_pot_01y2.o mod_baseYfunciones_01y2.o mod_photoini_01y2.o\
       mod_absorcion_01y2.o mod_lanczos_01y2.o\
       $dir3/main_boundlanz.f $dir/liboctdyn.f $dir3/liboctdynlanz.f \
       $dir/fit_general.f $dir/dipele_general.f $dir/coupling_general.f \
 -lfftw3

######################  distri.out & distriREAC.out
rm distri.out distriREAC.out
gfortran  -g -C -fcheck=all -o distri.out $dir2/distriwvp.f $dir/liboctdyn.f
gfortran  -g -C -fcheck=all -o distriREAC.out $dir2/distriREACwvp.f $dir/liboctdyn.f

######################  crp.out & cip.out & cipave.out
rm cip.out crp.out cipave.out
gfortran  -g -C -fcheck=all -o crp.out $dir2/CRP-fast.f $dir/liboctdyn.f
gfortran  -g -C -fcheck=all -o cip.out $dir2/CIP-fast.f $dir/liboctdyn.f

######################  rate.out & rates2s.out
rm rate.out rates2s.out
gfortran  -g -C -fcheck=all -o rate.out $dir2/rateFromSigma.f $dir/liboctdyn.f
gfortran  -g -C -fcheck=all -o rates2s.out $dir2/rate-s2s-fromCRP-extrapolation.f
gfortran  -g -C -fcheck=all -o inelastic-rates2s.out $dir2/rate-s2s-fromCIP-extrapolation.f

######################  sigma.out 
rm sigma.out
mpif77  -g -C -fcheck=all -o sigma.out $dir2/sigmaFromS2prod.f $dir/liboctdyn.f
######################  cheby-spectra.out 
rm cheby-spectra.out
gfortran  -g -C -fcheck=all -o cheby-spectra.out $dir4/cheby-spectra.f 

######################  Einstein.out
rm Einstein.out
mpif77 -o Einstein.out -g -C -fcheck=all  mod_gridYpara_01y2.o mod_baseYfunciones_01y2.o mod_pot_01y2.o mod_Hphi_01y2.o mod_photoini_01y2.o $photo/Acoefficient-bound-bound.f $dir/liboctdyn.f $dir/fit_general.f $dir/dipele_general.f $dir/coupling_general.f -lfftw3


echo "12 executable codes: mad3.out bndgrid.out Einstein.out distri.out  distriREAC.out  crp.out  cip.out inelastic-rates2s.out rate.out rates2s.out sigma.out  cheby-spectra.out"

echo "removing *.o and *.mod"
rm -f *.o *.mod
