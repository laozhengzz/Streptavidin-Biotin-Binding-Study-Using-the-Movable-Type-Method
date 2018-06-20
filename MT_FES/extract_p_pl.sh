#### Extract APo protein and PL free eneergy for each complex ######
#### $1 = RXn coord. X, $2 = Rxn coord. Y, $3 = Protein Torsion, $4 = Protein Intra non-bond., $5 = Protein solvation, $9 = PL Inter non-bond., $10 = PL solvation free energy
### Protein free energy = ($3 + ($4/60) + $5)
### Protein-Ligand free energy = ($3 + ($4/60) + ($9/6) + $10)


awk '{print $1, $2, ($3 + ($4/60) + $5), ($3 + ($4/60) + ($9/6) + $10)}' ../rn_coord_alldata_docked.txt > rncoord_p_pl.txt 
