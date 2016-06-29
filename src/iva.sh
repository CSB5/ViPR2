#!/bin/bash

# put into separate script because needs change of conda env which messes things up in snakefile (I think)

FQ1=$1
FQ2=$2
OUTD=$3
THREADS=$4
test -z "$FQ1" && exit 1
test -z "$FQ2" && exit 1
test -z "$OUTD" && exit 1
test -z "$THREADS" && exit 1

test -d $OUTD && rm -rf $OUTD
test -e "$FQ1" || exit 1
test -e "$FQ2" || exit 1

# WARN: we've seen weird cases where nested envs and the path meddling here
# completely messes things up


export PATH=/mnt/software/unstowable/anaconda/bin:$PATH
#echo "DEBUG: path before source activate iva: $PATH" 1>&2
source activate iva || echo "Couldn't source IVA" 1>&2

# source deactivate assumes PATH was prepended. the following will mess with the assumption

export PATH=/mnt/software/stow/KMC-2.2/bin/:$PATH
export PATH=/mnt/software/stow/smalt-0.7.6/bin/:$PATH
export PATH=/mnt/software/stow/samtools-0.1.19/bin/:$PATH
export PERL5LIB=/mnt/software/lib:"/mnt/software/lib/perl5/`perl -e 'printf "%vd", $^V;'`":/mnt/software/share/perl5:$PERL5LLIB
export PATH=/mnt/software/unstowable/mummer-3.23:$PATH

iva -f $FQ1 -r $FQ2 -t $THREADS $OUTD
echo "DEBUG: path before source activate iva: $PATH" 1>&2
