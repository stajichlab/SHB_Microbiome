#!/usr/bin/bash -l
#SBATCH -p batch -c 32 --mem 64gb -N 1 -n 1 --out logs/03_amptk_dada2.%A.log

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi

module load amptk
OUT=results/amptk.dada2
mkdir -p $OUT
IN=$(realpath input)
BASE=FD_202307_3Pool

pushd $OUT
if [ ! -f $BASE.demux.fq.gz ]; then
	echo "did not find $BASE.demux.fq.gz; rerun step pipeline/02_amptk_preprocess.sh"
	exit
fi
if [ ! -f $BASE.ASVs.fa ]; then
	amptk dada2 -i $BASE.demux.fq.gz --platform illumina -o $BASE --cpus $CPU --uchime_ref ITS
fi

if [ ! -f $BASE.filtered.otus.fa ]; then
  amptk filter -i $BASE.otu_table.txt -f $BASE.cluster.otus.fa -p 0.005
fi

if [ ! -f $BASE.otu_table.taxonomy.txt ]; then
  amptk taxonomy -f $BASE.ASVs.fa -i $BASE.final.txt -d ITS
fi


if [ ! -f $BASE.ASVs.guilds.txt ]; then
 amptk funguild -i $BASE.ASVs.otu_table.taxonomy.txt -o $BASE.ASVs.guilds.txt
fi


if [ ! -f $BASE.taxonomy.fix.txt ]; then
 perl ../../scripts/rdp_taxonmy2mat.pl<$BASE.taxonomy.txt>$BASE.taxonomy.fix.txt
fi

popd
