dir=~/Dropbox/MadWave3/MadWave3-v6/
pot=~/PES/H3+/crisSinglets

cp $dir/main_potini.f main_potini.f
cp $dir/mod_gridYpara_01y2.f mod_gridYpara_01y2.f
cp $dir/mod_baseYfunciones_01y2.f mod_baseYfunciones_01y2.f
cp $dir/mod_pot_01y2.f mod_pot_01y2.f
cp $dir/mod_photoini_01y2.f mod_photoini_01y2.f
cp $dir/mod_colini_01y2.f .
cp $dir/mod_Hphi_01y2.f .
cp $dir/liboctdyn.f .


cp $pot/pes-mad2.f  fit.f
cat $pot/diat-triat-withder-snglt.f >> fit.f
cp $pot/dim-H3psnglt.f90 .
cp $pot/largo-v2-withder.f90 .
cp $dir/dipele_general.f dipele.f

mpif77 -c dim-H3psnglt.f90 largo-v2-withder.f90

mpif77  -O3 -o pot.out mod_gridYpara_01y2.f mod_pot_01y2.f mod_baseYfunciones_01y2.f \
 mod_Hphi_01y2.f mod_photoini_01y2.f dim-H3psnglt.o largo-v2-withder.o\
 main_potini.f liboctdyn.f fit.f dipele.f \
 -lfftw3


