#!/bin/bash
set -euxo pipefail

# must load pg conda env first with the perigrine env 

fofn=$1 # must be a local path 
asm_dir=$2 # must be "./something"
threads=24

DATE=$(date +'%Y%m%d')
USER_N=$(whoami)
tmp_dir=/tmp/${USER_N}_${DATE}/pg_asms

# clean out the asm dir
#rm -rf $tmp_dir/$asm_dir

# make tmp dir and fofn there
mkdir -p $tmp_dir
cp $fofn  $tmp_dir

# move to tmp dir
pushd $tmp_dir
fofn=$(basename $fofn)

yes yes | pg_run.py asm  \
        $fofn $threads $threads $threads $threads $threads $threads $threads $threads $threads \
        --with-consensus --shimmer-r 3 --best_n_ovlp 8 \
        --output $asm_dir || echo "for some reason pg throws an error at the end of the piepline, ignore and check for the assembly"

if [ -f $asm_dir/p_ctg_cns.fa ]; then

        samtools faidx $asm_dir/p_ctg_cns.fa

        popd

        # remove a lot of stuff
		rm $fofn
        rm -rf $tmp_dir/$asm_dir/0-seqdb  $tmp_dir/$asm_dir/1-index  $tmp_dir/$asm_dir/2-ovlp

        # move results off of the SSD
        mv $tmp_dir/$asm_dir .

else
        echo "PG FAILED"
        popd
fi



