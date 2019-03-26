#!/bin/bash

for n in {1..100}
do
i=$(printf "%03d" $n)
cp fam_bicel.gro temp.gro 
cp /bio/rabhaska/new_bicelle/OUTPUT/Processed_fam_bicels/fam_bicel_w_nj_c_mol_fit_$i.xtc temp.xtc

python2 Fit_sph_fam_bicel.py
mv data.dat /bio/rabhaska/new_bicelle/OUTPUT/Processed_fam_bicels/fam_bicel_curv_$i.sdat

done

