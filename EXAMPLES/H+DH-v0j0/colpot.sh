dir=../../SRC
pot=../../PES/H3

mpif77  -O3 -o pot.out $dir/mod_gridYpara_01y2.f $dir/mod_pot_01y2.f $dir/mod_baseYfunciones_01y2.f \
 $dir/mod_Hphi_01y2.f $dir/mod_photoini_01y2.f \
 $dir/main_potini.f $dir/liboctdyn.f $pot/H3pot.f $pot/pothead.f $dir/dipele_general.f \
 $dir/coupling_general.f\
 -lfftw3

