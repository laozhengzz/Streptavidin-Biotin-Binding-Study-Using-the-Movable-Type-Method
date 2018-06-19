######################################################
###### run MATLAB in parallel ########################
##### written by Nupur Bansal ########################
#### 	Dated: 10/17/2016     ########################
######################################################

# Read in your protein PDB (remove hetatms)
# Extract the pocket part and save it as xxx_pocket.pdb
# This one is without ligand, but ligand can be added using same principles


## Run the MATLAB function to save all the information
## relevant to the running of loop generation code inc-
## -luding no. of residues, disuphide bonds, h-bonds,
## name of residues.etc All the info will be saved in
## chk.mat file (VERY IMPORTANT)

module load MATLAB-compiler/R2016a
#./mtflex_parallel_read "/mnt/scratch/bansalnu/MTflexb/parallel_10_17_2016/PDBs/1mk5_protein.pdb" "/mnt/scratch/bansalnu/MTflexb/parallel_10_17_2016/PDBs/1mk5_pocket.pdb" > out.txt;

#a=$(awk 'NR==1{print}' out.txt);  #Number of residues in pocket area
#echo $a

#rm out.txt

## Run for residue i to n/2, where n-total number of residues in a loop
## For each residue, use 30 predefined backbone sets possible.
## So, for 1st residue, there will be 30 jobs, with each backbone conf
## running parallely.Alongwith backbone, sidechains will also be generated.
## For second residue, maximum 30*30= 900 jobs (mostly will be less).
## But now, each individual job will have different backbone conf, while side-
## chains will be compacted in one run.

mm=12; # number of columns in tor.mat
a=4; ## for 1mk5/1mk5 (as n=8)

#############################################################################
######### FOR i = 1 
i=1;
for (( j=1; j<=$mm; j+=1 ));
	do
	mkdir ${i}_${j}/
	cd ${i}_${j}/
	cp ../*.m .
	cp ../tor.mat .
	cp ../chk.mat .
	cp ../res.pbs .
	cp ../mtflex_opposite_parallel .
	sed s:yyy:${i}:g res.pbs > run.pbs
	sed s:uuu:${j}:g run.pbs > run1.pbs
	qsub run1.pbs
	cd ../
done

### Seeing how many of the jobs are still running.
## Waiting for all of them to finish

ch=$(qstat -u bansalnu|grep 1mk5_opp_$i |wc -l);
while [ $ch -ge 1 ]
	do
       	sleep 1s
   	ch=$(qstat -u bansalnu|grep 1mk5_opp_$i |wc -l);
      	echo $ch
done
## Taking out all the chckpoint files generated and renaming them
d=1;
for (( j=1; j<=$mm; j+=1 ));
	do
	aa=$(awk 'NR==1{print}' ${i}_${j}/out.txt)
	if [ ! -z "$aa" ];
      		then
        	mv ${i}_${j}/chkpoint.mat chkpoint_${i}_${d}.mat
                d=$((d+1));
 	fi
	rm -r ${i}_${j}/
done
############## FOR i=1 end here ###################################


######## For i=2 to n/2   #########################################
for (( i=2; i<=$a; i+=1 ));
	do
	n1=$((i-1));
        echo $n1;
        Num1=$(ls chkpoint_${n1}_*.mat|wc -l)
        echo $Num1;
        for (( j=1; j<=$Num1; j+=1 ));
                do
		for (( k=1; k<=$mm; k+=1 ));
			do
			mkdir ${i}_${j}_${k}
			cd ${i}_${j}_${k}/
			cp ../*.m .
			cp ../tor.mat .
			cp ../chk.mat .
			cp ../chkpoint_${n1}_${j}.mat chkpoint.mat
                        cp ../res.pbs .
			cp ../mtflex_opposite_parallel .
			sed s:yyy:${i}:g res.pbs > run.pbs
                        sed s:uuu:${k}:g run.pbs > run1.pbs
			sed s:xxx:${j}:g run1.pbs > run2.pbs
			ch=$(qstat -u bansalnu |wc -l);
                        while [ $ch -ge 900 ]
                                do sleep 1s
                                ch=$(qstat -u bansalnu |wc -l);
                        done
                                qsub run2.pbs
                        cd ../
		done	
	done

	ch=$(qstat -u bansalnu|grep 1mk5_opp_$i |wc -l);
        while [ $ch -ge 1 ]
                do
                 sleep 1s
                ch=$(qstat -u bansalnu|grep 1mk5_opp_$i |wc -l);
                echo $ch
        done
        d=1;

	for (( j=1; j<=$Num1; j+=1 ));
                do
                for (( k=1; k<=$mm; k+=1 ));
                        do
                         aa=$(awk 'NR==1{print}' ${i}_${j}_${k}/out.txt)
                        if [ ! -z "$aa" ];
                                then
                                mv ${i}_${j}_${k}/chkpoint.mat chkpoint_${i}_${d}.mat
                                d=$((d+1));
                        fi
                        rm -r ${i}_${j}_${k}/
                done
        done
        rm chkpoint_${n1}_*.mat 
done
