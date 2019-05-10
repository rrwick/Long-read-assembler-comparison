#!/usr/bin/env bash

cd /projects/js66/individuals/ryan/Badread/assemblies

./completed_contigs.py > real_set_results
for s in SAMN10819801 SAMN10819805 SAMN10819807 SAMN10819813 SAMN10819815 SAMN10819847; do
    ref=../reads/real/"$s"/reference.fasta
    for p in nanopore pacbio; do
        for x in {00..09}; do
            reads=../reads/real/"$s"_"$p"_"$x".fastq.gz
            for assembly in real/"$s"_"$p"_"$x"_*.fasta.gz; do
                ls -l $assembly
                ./completed_contigs.py $assembly $reads $ref >> real_set_results
            done
        done
    done
done
