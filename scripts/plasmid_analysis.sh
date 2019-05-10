#!/usr/bin/env bash

cd /projects/js66/individuals/ryan/Badread/assemblies

./completed_contigs.py > plasmid_results
for i in {001..222}; do
    reads=../reads/plasmids/"$i".fastq.gz
    ref=../reads/plasmid_references/"$i".fasta
    for assembly in plasmids/"$i"_*.fasta.gz; do
        ./completed_contigs.py $assembly $reads $ref >> plasmid_results
    done
done
