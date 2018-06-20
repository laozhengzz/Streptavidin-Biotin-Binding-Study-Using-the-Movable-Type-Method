rm lowest_RMSDs.txt

sort -g -k2 rmsd_all_1mk5_back.dat > sorted_rmsd_all_1mk5_back.dat

sort -g -k2 rmsd_all_1swc_back.dat > sorted_rmsd_all_1swc_back.dat

awk 'FNR==2{print "RMSD_wrt_to_crystal_closed", $2}' sorted_rmsd_all_1mk5_back.dat > lowest_RMSDs.txt

awk 'FNR==2{print "RMSD_wrt_to_crystal_open", $2}' sorted_rmsd_all_1swc_back.dat >> lowest_RMSDs.txt

rm sorted_rmsd_all_1mk5_back.dat sorted_rmsd_all_1swc_back.dat
