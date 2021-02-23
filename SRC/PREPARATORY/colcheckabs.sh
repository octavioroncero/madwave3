dir=~/Dropbox/MadWave3/MadWave3-v6/PREPARATORY
pot=~/PES/H3+/crisSinglets

cp $pot/pes-mad2.f  fit.f

mpif77 -c $pot/dim-H3psnglt.f90 $pot/largo-v2-withder.f90

rm check.out

mpif77  -O3 -g -C -Wall -fcheck=all  -o check.out $dir/abscheck-in-r2.f\
  $dir/../liboctdyn.f \
  dim-H3psnglt.o largo-v2-withder.o fit.f \
 -lfftw3


