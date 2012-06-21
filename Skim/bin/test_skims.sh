codedir=$CMSSW_BASE/src/MitPhysics/Skim
com="root -b -l -q $codedir/macros/rootlogon.C"

export MIT_PROD_JSON=-
export MIT_PROD_OVERLAP=-1

run_skims() {
    . $codedir/config/$skim.env
    outdir=/scratch/dkralph/skims/$skim
    mkdir -p $outdir
    cd $outdir
    if echo "$outdir" | grep "/scratch/dkralph/" &>/dev/null; then rm -v $outdir/*{.root,.txt}; fi
    macro=$codedir/macros/$MIT_PROD_MACRO
    $com $macro+
    while read line; do
	if echo $line | grep '^#' &>/dev/null; then continue; fi
	dset=`echo $line | awk '{print $2}'`
	book=`echo $line | awk '{print $1}'`
	catalogDir=/home/cmsprod/catalog
	export MIT_PROD_JSON=`echo $line | awk '{print $8}'`
	$com $macro+\(\"0000\",\"noskim\",\"$dset\",\"$book\",\"$catalogDir\",\"$skim-$dset\",10000\) &> $outdir/$dset.txt &
    done < $codedir/config/$MIT_PROD_CFG.txt
    cd $codedir
}

count_skim() {
    outfile=~/public_html/skims/$skim.txt
    outdir=/scratch/dkralph/skims/$skim
    ntot_all=0
    nsel_all=0
    echo "${skim}:" | tee $outfile
    for txtfile in `ls $outdir/*.txt`; do
	nsel=`grep 'selected events' $txtfile | awk '{print $5'}`
	ntot=`grep 'total events' $txtfile | awk '{print $5'}`
	dset=`basename $txtfile | sed 's/.txt//'`
	printf "%20s%10d/%10d%10f\n" $dset $nsel $ntot `echo $nsel/$ntot | bc -l` | tee -a $outfile
	ntot_all=$(($ntot_all+$ntot))
	nsel_all=$(($nsel_all+$nsel))
    done
    printf "\n%20s%10d/%10d%10f\n" "total" $nsel_all $ntot_all `echo $nsel_all/$ntot_all | bc -l` | tee -a $outfile
    echo "-----------------------------------"
}

#for skim in h4l; do
#for skim in h4lzplusfake; do
#for skim in h4leletagprobe; do
#for skim in h4llf; do
#for skim in h4lmutagprobe; do
for skim in h4l h4lzplusfake h4llf h4leletagprobe h4lmutagprobe; do
#    run_skims &
    count_skim
done
