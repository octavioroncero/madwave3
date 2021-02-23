dir=~/SRC/XBC/MADWAVE3/MadWave3-v6/
dir=~/Dropbox/MadWave3/MadWave3-v6/
dir2=~/Dropbox/MadWave3/MadWave3-v6/AUXILIARy-CODES

pot=~/Dropbox/MadWave3/MadWave3-v6/

######################  mad3.out
rm mad3.out

      mpif77 -c  $dir/mod_gridYpara_01y2.f
      mpif77 -c  $dir/mod_pot_01y2.f
      mpif77 -c  $dir/mod_baseYfunciones_01y2.f
      mpif77 -c  $dir/mod_Hphi_01y2.f
      mpif77 -c  $dir/mod_photoini_01y2.f
      mpif77 -c  $dir/mod_colini_01y2.f
      mpif77 -c  $dir/mod_absorcion_01y2.f
      mpif77 -c  $dir/mod_flux_01y2.f
      mpif77 -c  $dir/mod_coortrans01_02.f

mpif77 -O3 -o mad3.out  mod_colini_01y2.o  mod_Hphi_01y2.o mod_gridYpara_01y2.o\
       mod_pot_01y2.o mod_baseYfunciones_01y2.o mod_photoini_01y2.o \
       mod_absorcion_01y2.o mod_flux_01y2.o mod_coortrans01_02.o\
       $dir/main_madwave3.f $dir/liboctdyn.f \
       $dir/fit_general.f $dir/dipele_general.f \
 -lfftw3

######################  distri.out & distriREAC.out
rm distri.out distriREAC.out
gfortran  -O3 -o distri.out $dir2/distriwvp.f $dir/liboctdyn.f
gfortran  -O3 -o distriREAC.out $dir2/distriREACwvp.f $dir/liboctdyn.f

######################  crp.out & cip.out
rm cip.out crp.out
gfortran  -O3 -o crp.out $dir2/CRP-fast.f $dir/liboctdyn.f
gfortran  -O3 -o cip.out $dir2/CIP-fast.f $dir/liboctdyn.f

######################  rate.out & rates2s.out
rm rate.out rates2s.out
gfortran  -O3 -o rate.out $dir2/rateFromSigma.f $dir/liboctdyn.f
gfortran  -O3 -o rates2s.out $dir2/rate-s2s-fromCRP-extrapolation.f

######################  sigma.out 
rm sigma.out

mpif77  -O3 -o sigma.out $dir2/sigmaFromS2prod.f $dir/liboctdyn.f

echo "programas generados 8: mad3.out distri.out  distriREAC.out  crp.out  cip.out rate.out rates2s.out sigma.out"
