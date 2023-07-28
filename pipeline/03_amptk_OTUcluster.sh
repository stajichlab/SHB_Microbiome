#!/usr/bin/bash -l
#SBATCH -p short -c 64 --mem 64gb -N 1 -n 1 --out logs/03_amptk_OTUcluster.%A.log

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
	echo "did not find $BASE.demux.fq.gz; rerun step pipeline/02_amptk_preprocess.sh"
	exit
fi

if [ ! -f  $BASE.cluster.otus.fa ]; then
	#amptk cluster -i $BASE.demux.fq.gz -o $BASE --uchime_ref ITS --usearch usearch9 --map_filtered -e 0.9 --cpus $CPU --unoise

	amptk cluster -i $BASE.demux.fq.gz -o $BASE --usearch usearch9 --map_filtered -e 0.9 --cpus $CPU
fi

if [ ! -f $BASE.filtered.otus.fa ]; then
 amptk filter -i $BASE.otu_table.txt -f $BASE.cluster.otus.fa -p 0.005
fi

if [ ! -f $BASE.otu_table.taxonomy.txt ]; then
 amptk taxonomy -f $BASE.filtered.otus.fa -i $BASE.final.txt -d ITS
fi

if [ ! -f $BASE.guilds.txt ]; then
 amptk funguild -i $BASE.otu_table.taxonomy.txt -o $BASE.guilds.txt
fi

if [ ! -f $BASE.taxonomy.fix.txt ]; then
 perl ../../scripts/rdp_taxonmy2mat.pl $BASE.taxonomy.txt > $BASE.taxonomy.fix.txt
fi

popd
