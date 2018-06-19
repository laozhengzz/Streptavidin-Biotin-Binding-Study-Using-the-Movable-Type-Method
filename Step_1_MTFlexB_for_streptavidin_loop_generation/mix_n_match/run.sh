for m in `seq 1 4126` ##4126`
do
  mkdir ${m}	
  cd ${m}/
  cp ../mix_two_halves .
  cp ../*.m .
  cp ../chk.mat .
  cp ../run.pbs .
  sed s:yyy:${m}:g run.pbs > run1.pbs
  ch=$(qstat -u bansalnu| wc -l);
  while [ $ch -ge 800 ]
  	do sleep 1s
    	ch=$(qstat -u bansalnu| wc -l);
  done
  qsub run1.pbs
  cd ../
done
