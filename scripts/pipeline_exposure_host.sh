#!/usr/bin/env bash
# Exposure, host-read and prebuilt-database analyses. Idempotent; outputs skipped if present.
# Stages: host (hamster genomes + minimap2 index), t1l (MRV type 1 Lang simulation + classification),
#         pconf (panel at confidence 0.1), oconf (14 runs at confidence 0.1 against the 17 prebuilt databases)
set -uo pipefail
K=/home/hm/kdcr
R3="$(dirname "$(readlink -f "$0")")"
W=$K/work_exposure; mkdir -p $W/logs
E=/home/hm/miniforge3/envs/kdcr/bin
KR=$E/kraken2
export PATH=$E:$PATH
MSL="msl29 msl30 msl31 msl32 msl33 msl34 msl35 msl36 msl37 msl38 msl39 msl40 msl41"
OFFICIAL="20201202 20210517 20220607 20220908 20221209 20230314 20230605 20231009 20240112 20240605 20240904 20241228 20250402 20250714 20251015 20260226 20260626"
RUNS14="SRR26352207 SRR26352206 SRR26352204 SRR26352215 SRR26352205 SRR26352203 SRR26352202 SRR26352217 SRR26352216 SRR26352195 SRR26352184 SRR26352212 SRR26352209 SRR26352208"
FQ=$K/data/fastq
R2W=$K/work_designs
classify() {  # db out_prefix r1 r2 [extra args]
  local db=$1 out=$2 r1=$3 r2=$4; shift 4
  [ -s ${out}_report.txt ] && return 0
  $KR --db $db --memory-mapping --paired --threads 6 "$@" --report ${out}_report.txt --output ${out}_output.txt $r1 $r2 > /dev/null 2>&1
}
dbpath() {  # design version
  case $1 in fixed) echo $K/kraken2_dbs_fixed/$2;; d1) echo $K/kraken2_dbs/$2;; d3) echo $R2W/design3/db_$2;;
             official) echo $R2W/official/k2_viral_$2;; current) echo $K/kraken2_dbs/current;; esac; }

host() {
  cd $W/host
  for z in picr chok1; do
    [ -s $z.fna ] && continue
    unzip -o -q $z.zip -d $z && cat $z/ncbi_dataset/data/*/*.fna > $z.fna && rm -rf $z
  done
  for z in picr chok1; do [ -s $z.mmi ] || minimap2 -x sr -t 6 -d $z.mmi $z.fna > $W/logs/mmi_$z.log 2>&1; done
}

t1l() {
  T=$W/t1l; mkdir -p $T/cls
  [ -s $T/t1l_2.fq ] || wgsim -N 5000 -1 150 -2 150 -e 0.01 -r 0 -R 0 -d 300 -s 30 -S 215 $R3/../sequences/mrv_T1L_10seg.fna $T/t1l_1.fq $T/t1l_2.fq > /dev/null 2>&1
  for d in fixed d1 d3; do for m in $MSL; do classify $(dbpath $d $m) $T/cls/${d}_$m $T/t1l_1.fq $T/t1l_2.fq; done; done
  for o in $OFFICIAL; do classify $(dbpath official $o) $T/cls/official_$o $T/t1l_1.fq $T/t1l_2.fq; done
  classify $(dbpath current x) $T/cls/current $T/t1l_1.fq $T/t1l_2.fq
  classify $(dbpath d3 msl41) $T/cls/d3_msl41_c0.1 $T/t1l_1.fq $T/t1l_2.fq --confidence 0.1
}

pconf() {
  P=$W/pconf; mkdir -p $P
  for m in $MSL; do classify $(dbpath d3 $m) $P/d3_$m $R2W/panel/panel_1.fq $R2W/panel/panel_2.fq --confidence 0.1 --report-zero-counts; done
  for o in $OFFICIAL; do classify $(dbpath official $o) $P/official_$o $R2W/panel/panel_1.fq $R2W/panel/panel_2.fq --confidence 0.1 --report-zero-counts; done
  classify $(dbpath current x) $P/current $R2W/panel/panel_1.fq $R2W/panel/panel_2.fq --confidence 0.1 --report-zero-counts
}

oconf() {
  O=$W/oconf; mkdir -p $O/sub
  for r in $RUNS14; do   # union of the reads classified at confidence 0 by any prebuilt database
    [ -s $O/sub/${r}_2.fq ] && continue
    for o in $OFFICIAL; do awk '$1=="C"{print $2}' $R2W/official/cls/$o/${r}_output.txt; done | sort -u > $O/sub/$r.ids
    for mate in 1 2; do
      awk 'NR==FNR{k[$1];next} FNR%4==1{split(substr($1,2),a,"/"); keep=(a[1] in k)} keep' $O/sub/$r.ids $FQ/${r}_${mate}.fastq > $O/sub/${r}_${mate}.fq
    done
  done
  for o in $OFFICIAL; do mkdir -p $O/$o
    for r in $RUNS14; do classify $(dbpath official $o) $O/$o/$r $O/sub/${r}_1.fq $O/sub/${r}_2.fq --confidence 0.1; done
  done
}

STAGES=${@:-host t1l pconf oconf}
for s in $STAGES; do echo "[$(date +%T)] start $s"; $s; echo "[$(date +%T)] end $s"; done
