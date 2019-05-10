#!/usr/bin/env bash

cd /projects/js66/individuals/ryan/Badread/assemblies

./assembly_stats.py > read_length_results
for i in 00250 00500 00750 01000 01250 01500 01750 02000 02250 02500 02750 03000 03250 03500 03750 04000 04250 04500 04750 05000 05250 05500 05750 06000 06250 06500 06750 07000 07250 07500 07750 08000 08250 08500 08750 09000 09250 09500 09750 10000 10250 10500 10750 11000 11250 11500 11750 12000 12250 12500 12750 13000 13250 13500 13750 14000 14250 14500 14750 15000 15250 15500 15750 16000 16250 16500 16750 17000 17250 17500 17750 18000 18250 18500 18750 19000 19250 19500 19750 20000 20250 20500 20750 21000 21250 21500 21750 22000 22250 22500 22750 23000 23250 23500 23750 24000 24250 24500 24750 25000 25250 25500 25750 26000 26250 26500 26750 27000 27250 27500 27750 28000 28250 28500 28750 29000 29250 29500 29750 30000 30250 30500 30750 31000 31250 31500 31750 32000 32250 32500 32750 33000 33250 33500 33750 34000 34250 34500 34750 35000 35250 35500 35750 36000 36250 36500 36750 37000 37250 37500 37750 38000 38250 38500 38750 39000 39250 39500 39750 40000 40250 40500 40750 41000 41250 41500 41750 42000 42250 42500 42750 43000 43250 43500 43750 44000 44250 44500 44750 45000 45250 45500 45750 46000 46250 46500 46750 47000 47250 47500 47750 48000 48250 48500 48750 49000 49250 49500 49750 50000; do
    r=../reads/read_length/"$i".fastq.gz
    echo $r
    for f in read_length/"$i"_*.fasta.gz; do
        >&2 echo "    "$f
        minimap2 -c -t 8 -r 10000 -g 10000 -x asm20 --eqx ref_tripled.fasta $f 2> /dev/null > "$f".paf
        ./assembly_stats.py $f $r "$f".paf 3565206 >> read_length_results
        rm "$f".paf
    done
done