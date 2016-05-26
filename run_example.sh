source /mnt/software/etc/gis.bashrc
#eval `/mnt/projects/rpd/init -d`
#use -a mnt-software
#use python-2.7
#source activate py3k || exit 1

outdir=$(mktemp -d --tmpdir=out) || exit 1
rmdir $outdir
echo "#Writing to $outdir"
FQ1=../vipr-test/MUX1051_PDH132_ATCACG_L002/PDH132_down1M_s1.fastq.gz
FQ2=../vipr-test/MUX1051_PDH132_ATCACG_L002/PDH132_down1M_s2.fastq.gz
REF=../vipr-test/ref/DENV_2.fa
PRIMERS=../vipr-test/D2_Primers.fa
cat<<EOF
./src/vipr2.py -1 $FQ1 -2 $FQ2 -o $outdir -r $REF --no-run -n PDH132 -p $PRIMERS
EOF