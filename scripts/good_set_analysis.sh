#!/usr/bin/env bash

cd /projects/js66/individuals/ryan/Badread/assemblies

./assembly_stats.py > good_set_results
for i in {000..499}; do
    r=../reads/good_sets/"$i".fastq.gz
    echo $r
    for f in good_sets/"$i"_*.fasta.gz; do
        >&2 echo "    "$f
        minimap2 -c -t 8 -r 10000 -g 10000 -x asm20 --eqx ref_tripled.fasta $f 2> /dev/null > "$f".paf
        ./assembly_stats.py $f $r "$f".paf 3565206 >> good_set_results
        rm "$f".paf
    done
done
