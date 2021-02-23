dir=~/SRC/XBC/MADWAVE3/MadWave3-v6/

cp $dir/main_madwave3.f .
cp $dir/mod_Hphi_mpiYopenmp.f mod_Hphi_01y2.f 
cp $dir/mod_gridYpara_01y2.f .
cp $dir/mod_baseYfunciones_01y2.f .
cp $dir/mod_pot_01y2.f .
cp $dir/mod_colini_01y2.f .
cp $dir/mod_absorcion_01y2.f .
cp $dir/mod_photoini_01y2.f .
cp $dir/mod_flux_01y2.f .
cp $dir/mod_coortrans01_02.f .
cp $dir/liboctdyn.f .

cp $dir/dipele_general.f dipele.f


cp ~/PES/H3/H3pot.f   fit.f


rm madYopenmp.out

module load fftw/3.3.6
module load gcc/6.3.0
module load openmpi


mpif77  -O3 -o madYopenmp.out -fopenmp mod_colini_01y2.f  mod_Hphi_01y2.f mod_gridYpara_01y2.f mod_pot_01y2.f mod_baseYfunciones_01y2.f mod_photoini_01y2.f \
 mod_absorcion_01y2.f mod_flux_01y2.f mod_coortrans01_02.f\
 main_madwave3.f liboctdyn.f fit.f dipele.f\
 -lfftw3


#mpif90.mpich -O3  dynamicmem.f madwave3.f liboctdyn.f fit.f  -lfftw3 -L/home/csic/fam/orv/fftw-3.1.2-icc/lib

