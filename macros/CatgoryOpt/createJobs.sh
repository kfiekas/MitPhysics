#! /bin/bash
#startCut=-0.5
startCut=-1
# this basestring nas to contain the number of categories, and the cuts applied before. There's no ordering needed, this is done in the macro
which=$1
baseString=$2
logdir=$3
mkdir -p $logdir
for i in `seq 401`
  do
  cutval=$(echo "scale=3; $startCut+($i-1)*0.005" | bc)
  echo "#! /bin/bash

/afs/cern.ch/user/m/mingyang/optimizeMVACats/runLimitAsym.sh "$which" "$baseString" "$cutval"

" > submit.cmd
  chmod +x submit.cmd
  
  bsub -q 1nh -o $logdir/job_${cutval}.log < ./submit.cmd

  rm submit.cmd
done
