sort -g -k2 rmsd_all_1mk5_back.dat > sorted_rmsd_all_1mk5_back.dat

awk '!seen[$2]++' sorted_rmsd_all_1mk5_back.dat > unique_backbone.txt
rm sorted_rmsd_all_1mk5_back.dat 
