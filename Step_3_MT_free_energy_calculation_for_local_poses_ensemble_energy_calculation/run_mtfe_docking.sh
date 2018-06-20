a=$(awk 'END{print NR}' /mnt/scratch/bansalnu/MT_conf_sort/uniq.txt);
echo $a
#for i in `seq 1 $a`
for i in `seq 3301 $a`
do
        m=$(awk 'FNR=='$i'{print $7}' /mnt/scratch/bansalnu/MT_conf_sort/uniq.txt);
        n1=$(awk 'FNR=='$i'{print $8}' /mnt/scratch/bansalnu/MT_conf_sort/uniq.txt);
	sed s:xxx:${m}:g run_mtfe_docking.pbs > run.pbs
	sed s:zzz:${i}:g run.pbs > run1.pbs
        sed s:yyy:${n1}:g run1.pbs > run_${i}.pbs
        ch=$(qstat -u bansalnu|wc -l);
        while [ $ch -ge 1000 ]
        do sleep 1s
                ch=$(qstat -u bansalnu|wc -l);
        done
        qsub run_${i}.pbs
done
