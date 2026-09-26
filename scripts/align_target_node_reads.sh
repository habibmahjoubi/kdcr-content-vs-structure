#!/usr/bin/env bash
# Alignment of every read counted at the RSV or MRV node in any database (reads.fa from exposure_host_prebuilt.py collect):
# megablast to the RefSeq viral library, blastn-short to the RSV and MRV genomes (longest exact match),
# minimap2 to the Chinese hamster genomes (CriGri-PICRH-1.0, GCF_003668045.3; CHO-K1 CriGri_1.0, GCF_000223135.1).
# minimap2 is run in 1-Gb index batches (-I 1G) to stay within memory.
set -uo pipefail
K=/home/hm/kdcr; S=$K/work_exposure/spur; H=$K/work_exposure/host
E=/home/hm/miniforge3/envs/kdcr/bin; export PATH=$E:$PATH
LIB=$K/kraken2_dbs/current/bulk/viral.1.1.genomic.taxid.fna
awk '/^>/{p=($0 ~ /kraken:taxid\|(11250|208893|208895|12814|538123|10886|351073)[ |]/ || $0 ~ /^>NC_0(38235|01803|01781)/)} p' $LIB > $S/targets.fna
grep -c ">" $S/targets.fna
blastn -task megablast -query $S/reads.fa -db $K/work/blastdb/viral_current -evalue 1e-5 -max_target_seqs 5 -max_hsps 1 \
       -outfmt '6 qseqid sseqid pident length evalue stitle' -num_threads 6 > $S/viral_megablast.tsv 2>/dev/null
blastn -task blastn-short -word_size 7 -evalue 1000 -query $S/reads.fa -subject $S/targets.fna -ungapped -perc_identity 100 \
       -outfmt '6 qseqid sseqid pident length' > $S/short_vs_targets.tsv 2>/dev/null
for g in picr chok1; do
  [ -s $S/host_$g.paf ] || minimap2 -x sr -c --secondary=no -I 1G -t 6 $H/$g.fna $S/reads.fa > $S/host_$g.paf 2> $S/host_$g.log
done
echo align_done
