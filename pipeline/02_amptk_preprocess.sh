#!/usr/bin/bash -l
#SBATCH -p short -c 64 --mem 64gb -N 1 -n 1 --out logs/02_amptk_preprocess.%A.log

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi

module load amptk
OUT=results/amptk
mkdir -p $OUT
IN=$(realpath input)
BASE=FD_202307_3Pool

pushd $OUT
if [ ! -f $BASE.demux.fq.gz ]; then
 	amptk illumina -i $IN --merge_method vsearch -f ITS1-F -r ITS2 --require_primer off \
		-o $BASE --usearch usearch10 --cpus $CPU --rescue_forward on --primer_mismatch 2 -l 250
fi

popd
