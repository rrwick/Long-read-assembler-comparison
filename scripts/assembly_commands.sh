# This file contains the actual assembly commands I ran. Unlike the README which shows assembly
# commands in an easy-to-read and generalisable format, this script is uglier (uses loops and
# variables) and is specific to the cluster I used (SLURM commands and directories). So it is NOT
# intended for others to run as-is.

# I gave assemblies up to 4 days of wall time. If they didn't finish after that, I deemed them failures.





# PREPPING DIRECTORIES
# -----------------------------------------------------------------------------
cd /projects/js66/individuals/ryan/Long_read_assembler_comparison
mkdir assemblies
cd assemblies
mkdir adapter_lengths chimeras glitches random_junk read_depth read_identity read_length good_sets plasmids real




# PROBLEM READ SETS
# -----------------------------------------------------------------------------
threads=4
for group in adapter_lengths chimeras glitches random_junk read_depth read_identity read_length good_sets; do
    cd /projects/js66/individuals/ryan/Long_read_assembler_comparison/assemblies/"$group"

    for r in /projects/js66/individuals/ryan/Long_read_assembler_comparison/reads/"$group"/*.fastq.gz; do
        temp_dir=/scratch/js66/ryan/badread_assemblies/"$RANDOM""$RANDOM"

        set_name=$(basename $r | sed s'|.fastq.gz||')

        # Ra 07364a1
        # There weren't any releases for Ra, so I'm using the commit checksum.
        working_dir="$temp_dir"_ra
        out_dir=$(pwd)
        assembly_commands="module load cmake; "
        assembly_commands+="date > "$set_name"_ra_07364a1.time; "
        assembly_commands+="mkdir "$working_dir"; "
        assembly_commands+="cd "$working_dir"; "
        assembly_commands+="ra -t "$threads" -x ont "$r" > "$out_dir"/"$set_name"_ra_07364a1.fasta; "
        assembly_commands+="cd "$out_dir"; "
        assembly_commands+="date >> "$set_name"_ra_07364a1.time; "
        assembly_commands+="gzip "$set_name"_ra_07364a1.fasta; "
        assembly_commands+="rm -r "$working_dir";"
        sbatch --job-name="$group"_"$set_name"_ra --wrap "$assembly_commands" --time=0-01:00:00 --partition=comp --account md52 --qos=normal --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_ra_07364a1.out --error "$set_name"_ra_07364a1.err
        sleep 0.2

        # Flye v2.4.2
        out_dir="$temp_dir"_flye
        temp_reads="$out_dir".fastq.gz  # Flye doesn't like commas in read filenames, so I make a copy to work on
        assembly_commands="module load python/2.7.15-gcc5; "
        assembly_commands+="cp "$r" "$temp_reads"; "
        assembly_commands+="date > "$set_name"_flye_v2.4.2.time; "
        assembly_commands+="/home/rwic0002/programs/Flye-2.4.2/bin/flye --nano-raw "$temp_reads" -g 3.565m -o "$out_dir" -t "$threads"; "  # I added the --plasmids option for the plasmid test sets
        assembly_commands+="date >> "$set_name"_flye_v2.4.2.time; "
        assembly_commands+="cp "$out_dir"/scaffolds.fasta "$set_name"_flye_v2.4.2.fasta; "
        assembly_commands+="cp "$out_dir"/assembly_graph.gfa "$set_name"_flye_v2.4.2.gfa; "
        assembly_commands+="gzip "$set_name"_flye_v2.4.2.fasta "$set_name"_flye_v2.4.2.gfa; "
        assembly_commands+="rm -r "$out_dir"; "
        assembly_commands+="rm "$temp_reads";"
        sbatch --job-name="$group"_"$set_name"_flye --wrap "$assembly_commands" --time=0-04:00:00 --partition=comp --account md52 --qos=normal --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_flye_v2.4.2.out --error "$set_name"_flye_v2.4.2.err
        sleep 0.2

        # Wtdbg2 v2.4
        # This one is a bit more complex to run, with assembly, consensus and polishing steps.
        out_dir="$temp_dir"_wtdbg2
        out_prefix="$out_dir"/wtdbg2
        assembly_commands="module load samtools/1.9; "
        assembly_commands+="date > "$set_name"_wtdbg2_v2.4.time; "
        assembly_commands+="mkdir "$out_dir"; "
        assembly_commands+="wtdbg2 -x nanopore -g 3.565m -i "$r" -t "$threads" -fo "$out_prefix"; "
        assembly_commands+="wtpoa-cns -t "$threads" -i "$out_prefix".ctg.lay.gz -fo "$out_prefix".ctg.fa; "
        assembly_commands+="minimap2 -t "$threads" -x map-ont -a "$out_prefix".ctg.fa "$r" | samtools sort > "$out_prefix".ctg.bam; "
        assembly_commands+="samtools view "$out_prefix".ctg.bam | wtpoa-cns -t "$threads" -d "$out_prefix".ctg.fa -i - -fo "$out_prefix".ctg.2nd.fa; "
        assembly_commands+="date >> "$set_name"_wtdbg2_v2.4.time; "
        assembly_commands+="cp "$out_prefix".ctg.2nd.fa "$set_name"_wtdbg2_v2.4.fasta; "
        assembly_commands+="gzip "$set_name"_wtdbg2_v2.4.fasta; "
        assembly_commands+="rm -r "$out_dir";"
        sbatch --job-name="$group"_"$set_name"_wtdbg2 --wrap "$assembly_commands" --time=0-00:30:00 --partition=comp --account md52 --qos=shortq --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_wtdbg2_v2.4.out --error "$set_name"_wtdbg2_v2.4.err
        sleep 0.2

        # Unicycler v0.4.7
        out_dir="$temp_dir"_unicycler
        assembly_commands="module purge; "
        assembly_commands+="module load python/3.6.6-gcc5; "
        assembly_commands+="module load blast/2.7.1; "
        assembly_commands+="date > "$set_name"_unicycler_v0.4.7.time; "
        assembly_commands+="/home/rwic0002/programs/Unicycler-0.4.7/unicycler-runner.py -l "$r" -o "$out_dir" -t "$threads"; "
        assembly_commands+="date >> "$set_name"_unicycler_v0.4.7.time; "
        assembly_commands+="cp "$out_dir"/assembly.fasta "$set_name"_unicycler_v0.4.7.fasta; "
        assembly_commands+="cp "$out_dir"/assembly.gfa "$set_name"_unicycler_v0.4.7.gfa; "
        assembly_commands+="gzip "$set_name"_unicycler_v0.4.7.fasta "$set_name"_unicycler_v0.4.7.gfa; "
        assembly_commands+="rm -r "$out_dir";"
        sbatch --job-name="$group"_"$set_name"_unicycler --wrap "$assembly_commands" --time=0-02:00:00 --partition=comp --account md52 --qos=normal --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_unicycler_v0.4.7.out --error "$set_name"_unicycler_v0.4.7.err
        sleep 0.2

        # Canu v1.8
        out_dir="$temp_dir"_canu
        assembly_commands="date > "$set_name"_canu_v1.8.time; "
        assembly_commands+="/home/rwic0002/programs/canu-1.8/Linux-amd64/bin/canu -p canu -d "$out_dir" genomeSize=3.565m stopOnLowCoverage=0 useGrid=false minThreads="$threads" maxThreads="$threads" maxMemory=31 -nanopore-raw "$r"; "
        assembly_commands+="date >> "$set_name"_canu_v1.8.time; "
        assembly_commands+="cp "$out_dir"/canu.contigs.fasta "$set_name"_canu_v1.8.fasta; "
        assembly_commands+="cp "$out_dir"/canu.contigs.gfa "$set_name"_canu_v1.8.gfa; "
        assembly_commands+="gzip "$set_name"_canu_v1.8.fasta "$set_name"_canu_v1.8.gfa; "
        assembly_commands+="rm -r "$out_dir";"
        sbatch --job-name="$group"_"$set_name"_canu --wrap "$assembly_commands" --time=0-18:00:00 --partition=comp --account md52 --qos=normal --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_canu_v1.8.out --error "$set_name"_canu_v1.8.err
        sleep 0.2

        # # SMARTdenovo 5cc1356
        # # There weren't any releases for SMARTdenovo, so I'm using the commit checksum.
        # out_dir="$temp_dir"_smartdenovo
        # out_prefix="$out_dir"/smartdenovo
        # assembly_commands="mkdir "$out_dir"; "
        # assembly_commands+="date > "$set_name"_smartdenovo_5cc1356.time; "
        # assembly_commands+="/home/rwic0002/programs/smartdenovo/smartdenovo.pl -p "$out_prefix" -c 1 "$r" -t "$threads" > "$out_prefix".mak; "
        # assembly_commands+="make -f "$out_prefix".mak; "
        # assembly_commands+="date >> "$set_name"_smartdenovo_5cc1356.time; "
        # assembly_commands+="cp "$out_prefix".cns "$set_name"_smartdenovo_5cc1356.fasta; "
        # assembly_commands+="gzip "$set_name"_smartdenovo_5cc1356.fasta; "
        # assembly_commands+="rm -r "$out_dir";"
        # sbatch --job-name="$group"_"$set_name"_smartdenovo --wrap "$assembly_commands" --time=0-04:00:00 --partition=comp --account md52 --qos=normal --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_smartdenovo_5cc1356.out --error "$set_name"_smartdenovo_5cc1356.err
        # sleep 0.2

    done
done






# PLASMID READ SETS
# -----------------------------------------------------------------------------
# The commands for launching plasmids are a bit different than the problem read sets:
# 1) Some of the assemblers need a genome size and unlike for the problem sets (which all used the
#    same one-chromosome reference), this can now differ. So I used fast_count to get the genome
#    size (plasmids included) from the reference.
# 2) I tried Flye in four configurations: with and without --plasmids and --meta in all
#    combinations.
threads=4
cd /projects/js66/individuals/ryan/Long_read_assembler_comparison/assemblies/plasmids

for x in {001..215}; do
    r=/projects/js66/individuals/ryan/Long_read_assembler_comparison/reads/plasmids/"$x".fastq.gz

    ref=/projects/js66/individuals/ryan/Long_read_assembler_comparison/reads/plasmid_references/"$x".fasta
    genome_size=$(echo "scale=3; "$(fast_count $ref | cut -f3)"/1000000" | bc -l)"m"

    temp_dir=/scratch/js66/ryan/badread_assemblies/"$RANDOM""$RANDOM"

    set_name=$(basename $r | sed s'|.fastq.gz||')

    # Ra 07364a1
    # There weren't any releases for Ra, so I'm using the commit checksum.
    working_dir="$temp_dir"_ra
    out_dir=$(pwd)
    assembly_commands="module load cmake; "
    assembly_commands+="date > "$set_name"_ra_07364a1.time; "
    assembly_commands+="mkdir "$working_dir"; "
    assembly_commands+="cd "$working_dir"; "
    assembly_commands+="ra -t "$threads" -x ont "$r" > "$out_dir"/"$set_name"_ra_07364a1.fasta; "
    assembly_commands+="cd "$out_dir"; "
    assembly_commands+="date >> "$set_name"_ra_07364a1.time; "
    assembly_commands+="gzip "$set_name"_ra_07364a1.fasta; "
    assembly_commands+="rm -r "$working_dir";"
    sbatch --job-name=plasmids_"$set_name"_ra --wrap "$assembly_commands" --time=0-01:00:00 --partition=comp --account md52 --qos=normal --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_ra_07364a1.out --error "$set_name"_ra_07364a1.err
    sleep 0.2

    # Flye v2.4.2
    out_dir="$temp_dir"_flye
    assembly_commands="module load python/2.7.15-gcc5; "
    assembly_commands+="date > "$set_name"_flye_v2.4.2.time; "
    assembly_commands+="/home/rwic0002/programs/Flye-2.4.2/bin/flye --nano-raw "$r" -g "$genome_size" -o "$out_dir" -t "$threads"; "  # I added the --plasmids option for the plasmid test sets
    assembly_commands+="date >> "$set_name"_flye_v2.4.2.time; "
    assembly_commands+="cp "$out_dir"/scaffolds.fasta "$set_name"_flye_v2.4.2.fasta; "
    assembly_commands+="cp "$out_dir"/assembly_graph.gfa "$set_name"_flye_v2.4.2.gfa; "
    assembly_commands+="gzip "$set_name"_flye_v2.4.2.fasta "$set_name"_flye_v2.4.2.gfa; "
    assembly_commands+="rm -r "$out_dir"; "
    sbatch --job-name=plasmids_"$set_name"_flye --wrap "$assembly_commands" --time=0-04:00:00 --partition=comp --account md52 --qos=normal --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_flye_v2.4.2.out --error "$set_name"_flye_v2.4.2.err
    sleep 0.2

    # Flye v2.4.2 --plasmids
    out_dir="$temp_dir"_flye_plasmids
    assembly_commands="module load python/2.7.15-gcc5; "
    assembly_commands+="date > "$set_name"_flye_v2.4.2_plasmids.time; "
    assembly_commands+="/home/rwic0002/programs/Flye-2.4.2/bin/flye --plasmids --nano-raw "$r" -g "$genome_size" -o "$out_dir" -t "$threads"; "
    assembly_commands+="date >> "$set_name"_flye_v2.4.2_plasmids.time; "
    assembly_commands+="cp "$out_dir"/scaffolds.fasta "$set_name"_flye_v2.4.2_plasmids.fasta; "
    assembly_commands+="cp "$out_dir"/assembly_graph.gfa "$set_name"_flye_v2.4.2_plasmids.gfa; "
    assembly_commands+="gzip "$set_name"_flye_v2.4.2_plasmids.fasta "$set_name"_flye_v2.4.2_plasmids.gfa; "
    assembly_commands+="rm -r "$out_dir"; "
    sbatch --job-name=plasmids_"$set_name"_flye_plasmids --wrap "$assembly_commands" --time=0-04:00:00 --partition=comp --account md52 --qos=normal --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_flye_v2.4.2_plasmids.out --error "$set_name"_flye_v2.4.2_plasmids.err
    sleep 0.2

    # Flye v2.4.2 --meta
    out_dir="$temp_dir"_flye_meta
    assembly_commands="module load python/2.7.15-gcc5; "
    assembly_commands+="date > "$set_name"_flye_v2.4.2_meta.time; "
    assembly_commands+="/home/rwic0002/programs/Flye-2.4.2/bin/flye --meta --nano-raw "$r" -g "$genome_size" -o "$out_dir" -t "$threads"; "
    assembly_commands+="date >> "$set_name"_flye_v2.4.2_meta.time; "
    assembly_commands+="cp "$out_dir"/scaffolds.fasta "$set_name"_flye_v2.4.2_meta.fasta; "
    assembly_commands+="cp "$out_dir"/assembly_graph.gfa "$set_name"_flye_v2.4.2_meta.gfa; "
    assembly_commands+="gzip "$set_name"_flye_v2.4.2_meta.fasta "$set_name"_flye_v2.4.2_meta.gfa; "
    assembly_commands+="rm -r "$out_dir"; "
    sbatch --job-name=plasmids_"$set_name"_flye_meta --wrap "$assembly_commands" --time=0-04:00:00 --partition=comp --account md52 --qos=normal --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_flye_v2.4.2_meta.out --error "$set_name"_flye_v2.4.2_meta.err
    sleep 0.2

    # Flye v2.4.2 --plasmids --meta
    out_dir="$temp_dir"_flye_plasmids_meta
    assembly_commands="module load python/2.7.15-gcc5; "
    assembly_commands+="date > "$set_name"_flye_v2.4.2_plasmids_meta.time; "
    assembly_commands+="/home/rwic0002/programs/Flye-2.4.2/bin/flye --plasmids --meta --nano-raw "$r" -g "$genome_size" -o "$out_dir" -t "$threads"; "
    assembly_commands+="date >> "$set_name"_flye_v2.4.2_plasmids_meta.time; "
    assembly_commands+="cp "$out_dir"/scaffolds.fasta "$set_name"_flye_v2.4.2_plasmids_meta.fasta; "
    assembly_commands+="cp "$out_dir"/assembly_graph.gfa "$set_name"_flye_v2.4.2_plasmids_meta.gfa; "
    assembly_commands+="gzip "$set_name"_flye_v2.4.2_plasmids_meta.fasta "$set_name"_flye_v2.4.2_plasmids_meta.gfa; "
    assembly_commands+="rm -r "$out_dir"; "
    sbatch --job-name=plasmids_"$set_name"_flye_plasmids_meta --wrap "$assembly_commands" --time=0-04:00:00 --partition=comp --account md52 --qos=normal --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_flye_v2.4.2_plasmids_meta.out --error "$set_name"_flye_v2.4.2_plasmids_meta.err
    sleep 0.2

    # Wtdbg2 v2.4
    # This one is a bit more complex to run, with assembly, consensus and polishing steps.
    out_dir="$temp_dir"_wtdbg2
    out_prefix="$out_dir"/wtdbg2
    assembly_commands="module load samtools/1.9; "
    assembly_commands+="date > "$set_name"_wtdbg2_v2.4.time; "
    assembly_commands+="mkdir "$out_dir"; "
    assembly_commands+="wtdbg2 -x nanopore -g "$genome_size" -i "$r" -t "$threads" -fo "$out_prefix"; "  # I added "-L 0" to the command for the read length tests
    assembly_commands+="wtpoa-cns -t "$threads" -i "$out_prefix".ctg.lay.gz -fo "$out_prefix".ctg.fa; "
    assembly_commands+="minimap2 -t "$threads" -x map-ont -a "$out_prefix".ctg.fa "$r" | samtools sort > "$out_prefix".ctg.bam; "
    assembly_commands+="samtools view "$out_prefix".ctg.bam | wtpoa-cns -t "$threads" -d "$out_prefix".ctg.fa -i - -fo "$out_prefix".ctg.2nd.fa; "
    assembly_commands+="date >> "$set_name"_wtdbg2_v2.4.time; "
    assembly_commands+="cp "$out_prefix".ctg.2nd.fa "$set_name"_wtdbg2_v2.4.fasta; "
    assembly_commands+="gzip "$set_name"_wtdbg2_v2.4.fasta; "
    assembly_commands+="rm -r "$out_dir";"
    sbatch --job-name=plasmids_"$set_name"_wtdbg2 --wrap "$assembly_commands" --time=0-00:30:00 --partition=comp --account md52 --qos=shortq --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_wtdbg2_v2.4.out --error "$set_name"_wtdbg2_v2.4.err
    sleep 0.2

    # Unicycler v0.4.7
    out_dir="$temp_dir"_unicycler
    assembly_commands="module purge; "
    assembly_commands+="module load python/3.6.6-gcc5; "
    assembly_commands+="module load blast/2.7.1; "
    assembly_commands+="date > "$set_name"_unicycler_v0.4.7.time; "
    assembly_commands+="/home/rwic0002/programs/Unicycler-0.4.7/unicycler-runner.py -l "$r" -o "$out_dir" -t "$threads"; "
    assembly_commands+="date >> "$set_name"_unicycler_v0.4.7.time; "
    assembly_commands+="cp "$out_dir"/assembly.fasta "$set_name"_unicycler_v0.4.7.fasta; "
    assembly_commands+="cp "$out_dir"/assembly.gfa "$set_name"_unicycler_v0.4.7.gfa; "
    assembly_commands+="gzip "$set_name"_unicycler_v0.4.7.fasta "$set_name"_unicycler_v0.4.7.gfa; "
    assembly_commands+="rm -r "$out_dir";"
    sbatch --job-name=plasmids_"$set_name"_unicycler --wrap "$assembly_commands" --time=0-02:00:00 --partition=comp --account md52 --qos=normal --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_unicycler_v0.4.7.out --error "$set_name"_unicycler_v0.4.7.err
    sleep 0.2

    # Canu v1.8
    out_dir="$temp_dir"_canu
    assembly_commands="date > "$set_name"_canu_v1.8.time; "
    assembly_commands+="/home/rwic0002/programs/canu-1.8/Linux-amd64/bin/canu -p canu -d "$out_dir" genomeSize="$genome_size" stopOnLowCoverage=0 useGrid=false minThreads="$threads" maxThreads="$threads" maxMemory=31 -nanopore-raw "$r"; "
    assembly_commands+="date >> "$set_name"_canu_v1.8.time; "
    assembly_commands+="cp "$out_dir"/canu.contigs.fasta "$set_name"_canu_v1.8.fasta; "
    assembly_commands+="cp "$out_dir"/canu.contigs.gfa "$set_name"_canu_v1.8.gfa; "
    assembly_commands+="gzip "$set_name"_canu_v1.8.fasta "$set_name"_canu_v1.8.gfa; "
    assembly_commands+="rm -r "$out_dir";"
    sbatch --job-name=plasmids_"$set_name"_canu --wrap "$assembly_commands" --time=4-00:00:00 --partition=comp --account md52 --qos=normal --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_canu_v1.8.out --error "$set_name"_canu_v1.8.err
    sleep 0.2
done






# REAL READ SETS
# -----------------------------------------------------------------------------
# Like was the case for the plasmid test sets, these commands need to get the correct reference
# genome size. Also, some of the commands need to change a bit based on whether it is a PacBio or
# a Nanopore read set being assembled.
threads=4
cd /projects/js66/individuals/ryan/Long_read_assembler_comparison/assemblies/real

for f in SAMN10819801 SAMN10819805 SAMN10819807 SAMN10819813 SAMN10819815 SAMN10819847; do

    ref=/projects/js66/individuals/ryan/Long_read_assembler_comparison/reads/real/"$f"/reference.fasta
    genome_size=$(echo "scale=3; "$(fast_count $ref | cut -f3)"/1000000" | bc -l)"m"

    for s in nanopore pacbio; do

        if [ $s = "pacbio" ]; then
           wtdbg2_preset="rsII"
           minimap2_preset="map-pb"
           canu_input="-pacbio-raw"
           ra_preset="pb"
           flye_input="--pacbio-raw"
        else
           wtdbg2_preset="nanopore"
           minimap2_preset="map-ont"
           canu_input="-nanopore-raw"
           ra_preset="ont"
           flye_input="--nano-raw"
        fi

        for x in {00..09}; do
            r=/projects/js66/individuals/ryan/Long_read_assembler_comparison/reads/real/"$f"_"$s"_"$x".fastq.gz
            temp_dir=/scratch/js66/ryan/badread_assemblies/"$RANDOM""$RANDOM"

            set_name=$(basename $r | sed s'|.fastq.gz||')

            # Ra 07364a1
            # There weren't any releases for Ra, so I'm using the commit checksum.
            working_dir="$temp_dir"_ra
            out_dir=$(pwd)
            assembly_commands="module load cmake; "
            assembly_commands+="date > "$set_name"_ra_07364a1.time; "
            assembly_commands+="mkdir "$working_dir"; "
            assembly_commands+="cd "$working_dir"; "
            assembly_commands+="ra -t "$threads" -x "$ra_preset" "$r" > "$out_dir"/"$set_name"_ra_07364a1.fasta; "
            assembly_commands+="cd "$out_dir"; "
            assembly_commands+="date >> "$set_name"_ra_07364a1.time; "
            assembly_commands+="gzip "$set_name"_ra_07364a1.fasta; "
            assembly_commands+="rm -r "$working_dir";"
            sbatch --job-name=real_"$set_name"_ra --wrap "$assembly_commands" --time=0-01:00:00 --partition=comp --account md52 --qos=normal --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_ra_07364a1.out --error "$set_name"_ra_07364a1.err
            sleep 0.2

            # Flye v2.4.2
            out_dir="$temp_dir"_flye
            assembly_commands="module load python/2.7.15-gcc5; "
            assembly_commands+="date > "$set_name"_flye_v2.4.2.time; "
            assembly_commands+="/home/rwic0002/programs/Flye-2.4.2/bin/flye "$flye_input" "$r" -g "$genome_size" -o "$out_dir" -t "$threads"; "  # I added the --plasmids option for the plasmid test sets
            assembly_commands+="date >> "$set_name"_flye_v2.4.2.time; "
            assembly_commands+="cp "$out_dir"/scaffolds.fasta "$set_name"_flye_v2.4.2.fasta; "
            assembly_commands+="cp "$out_dir"/assembly_graph.gfa "$set_name"_flye_v2.4.2.gfa; "
            assembly_commands+="gzip "$set_name"_flye_v2.4.2.fasta "$set_name"_flye_v2.4.2.gfa; "
            assembly_commands+="rm -r "$out_dir"; "
            sbatch --job-name=real_"$set_name"_flye --wrap "$assembly_commands" --time=0-08:00:00 --partition=comp --account md52 --qos=normal --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_flye_v2.4.2.out --error "$set_name"_flye_v2.4.2.err
            sleep 0.2

            # Flye v2.4.2 --plasmids
            out_dir="$temp_dir"_flye_plasmids
            assembly_commands="module load python/2.7.15-gcc5; "
            assembly_commands+="date > "$set_name"_flye_v2.4.2_plasmids.time; "
            assembly_commands+="/home/rwic0002/programs/Flye-2.4.2/bin/flye --plasmids "$flye_input" "$r" -g "$genome_size" -o "$out_dir" -t "$threads"; "
            assembly_commands+="date >> "$set_name"_flye_v2.4.2_plasmids.time; "
            assembly_commands+="cp "$out_dir"/scaffolds.fasta "$set_name"_flye_v2.4.2_plasmids.fasta; "
            assembly_commands+="cp "$out_dir"/assembly_graph.gfa "$set_name"_flye_v2.4.2_plasmids.gfa; "
            assembly_commands+="gzip "$set_name"_flye_v2.4.2_plasmids.fasta "$set_name"_flye_v2.4.2_plasmids.gfa; "
            assembly_commands+="rm -r "$out_dir"; "
            sbatch --job-name=real_"$set_name"_flye_plasmids --wrap "$assembly_commands" --time=0-08:00:00 --partition=comp --account md52 --qos=normal --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_flye_v2.4.2_plasmids.out --error "$set_name"_flye_v2.4.2_plasmids.err
            sleep 0.2

            # Flye v2.4.2 --meta
            out_dir="$temp_dir"_flye_meta
            assembly_commands="module load python/2.7.15-gcc5; "
            assembly_commands+="date > "$set_name"_flye_v2.4.2_meta.time; "
            assembly_commands+="/home/rwic0002/programs/Flye-2.4.2/bin/flye --meta "$flye_input" "$r" -g "$genome_size" -o "$out_dir" -t "$threads"; "
            assembly_commands+="date >> "$set_name"_flye_v2.4.2_meta.time; "
            assembly_commands+="cp "$out_dir"/scaffolds.fasta "$set_name"_flye_v2.4.2_meta.fasta; "
            assembly_commands+="cp "$out_dir"/assembly_graph.gfa "$set_name"_flye_v2.4.2_meta.gfa; "
            assembly_commands+="gzip "$set_name"_flye_v2.4.2_meta.fasta "$set_name"_flye_v2.4.2_meta.gfa; "
            assembly_commands+="rm -r "$out_dir"; "
            sbatch --job-name=real_"$set_name"_flye_meta --wrap "$assembly_commands" --time=0-08:00:00 --partition=comp --account md52 --qos=normal --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_flye_v2.4.2_meta.out --error "$set_name"_flye_v2.4.2_meta.err
            sleep 0.2

            # Flye v2.4.2 --plasmids --meta
            out_dir="$temp_dir"_flye_plasmids_meta
            assembly_commands="module load python/2.7.15-gcc5; "
            assembly_commands+="date > "$set_name"_flye_v2.4.2_plasmids_meta.time; "
            assembly_commands+="/home/rwic0002/programs/Flye-2.4.2/bin/flye --plasmids --meta "$flye_input" "$r" -g "$genome_size" -o "$out_dir" -t "$threads"; "
            assembly_commands+="date >> "$set_name"_flye_v2.4.2_plasmids_meta.time; "
            assembly_commands+="cp "$out_dir"/scaffolds.fasta "$set_name"_flye_v2.4.2_plasmids_meta.fasta; "
            assembly_commands+="cp "$out_dir"/assembly_graph.gfa "$set_name"_flye_v2.4.2_plasmids_meta.gfa; "
            assembly_commands+="gzip "$set_name"_flye_v2.4.2_plasmids_meta.fasta "$set_name"_flye_v2.4.2_plasmids_meta.gfa; "
            assembly_commands+="rm -r "$out_dir"; "
            sbatch --job-name=real_"$set_name"_flye_plasmids_meta --wrap "$assembly_commands" --time=0-08:00:00 --partition=comp --account md52 --qos=normal --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_flye_v2.4.2_plasmids_meta.out --error "$set_name"_flye_v2.4.2_plasmids_meta.err
            sleep 0.2

            # Wtdbg2 v2.4
            # This one is a bit more complex to run, with assembly, consensus and polishing steps.
            out_dir="$temp_dir"_wtdbg2
            out_prefix="$out_dir"/wtdbg2
            assembly_commands="module load samtools/1.9; "
            assembly_commands+="date > "$set_name"_wtdbg2_v2.4.time; "
            assembly_commands+="mkdir "$out_dir"; "
            assembly_commands+="wtdbg2 -x "$wtdbg2_preset" -g "$genome_size" -i "$r" -t "$threads" -fo "$out_prefix"; "
            assembly_commands+="wtpoa-cns -t "$threads" -i "$out_prefix".ctg.lay.gz -fo "$out_prefix".ctg.fa; "
            assembly_commands+="minimap2 -t "$threads" -x "$minimap2_preset" -a "$out_prefix".ctg.fa "$r" | samtools sort > "$out_prefix".ctg.bam; "
            assembly_commands+="samtools view "$out_prefix".ctg.bam | wtpoa-cns -t "$threads" -d "$out_prefix".ctg.fa -i - -fo "$out_prefix".ctg.2nd.fa; "
            assembly_commands+="date >> "$set_name"_wtdbg2_v2.4.time; "
            assembly_commands+="cp "$out_prefix".ctg.2nd.fa "$set_name"_wtdbg2_v2.4.fasta; "
            assembly_commands+="gzip "$set_name"_wtdbg2_v2.4.fasta; "
            assembly_commands+="rm -r "$out_dir";"
            sbatch --job-name=real_"$set_name"_wtdbg2 --wrap "$assembly_commands" --time=0-00:30:00 --partition=comp --account md52 --qos=shortq --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_wtdbg2_v2.4.out --error "$set_name"_wtdbg2_v2.4.err
            sleep 0.2

            # Unicycler v0.4.7
            out_dir="$temp_dir"_unicycler
            assembly_commands="module purge; "
            assembly_commands+="module load python/3.6.6-gcc5; "
            assembly_commands+="module load blast/2.7.1; "
            assembly_commands+="date > "$set_name"_unicycler_v0.4.7.time; "
            assembly_commands+="/home/rwic0002/programs/Unicycler-0.4.7/unicycler-runner.py -l "$r" -o "$out_dir" -t "$threads"; "
            assembly_commands+="date >> "$set_name"_unicycler_v0.4.7.time; "
            assembly_commands+="cp "$out_dir"/assembly.fasta "$set_name"_unicycler_v0.4.7.fasta; "
            assembly_commands+="cp "$out_dir"/assembly.gfa "$set_name"_unicycler_v0.4.7.gfa; "
            assembly_commands+="gzip "$set_name"_unicycler_v0.4.7.fasta "$set_name"_unicycler_v0.4.7.gfa; "
            assembly_commands+="rm -r "$out_dir";"
            sbatch --job-name=real_"$set_name"_unicycler --wrap "$assembly_commands" --time=0-02:00:00 --partition=comp --account md52 --qos=normal --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_unicycler_v0.4.7.out --error "$set_name"_unicycler_v0.4.7.err
            sleep 0.2

            # Canu v1.8
            out_dir="$temp_dir"_canu
            assembly_commands="date > "$set_name"_canu_v1.8.time; "
            assembly_commands+="/home/rwic0002/programs/canu-1.8/Linux-amd64/bin/canu -p canu -d "$out_dir" genomeSize="$genome_size" stopOnLowCoverage=0 useGrid=false minThreads="$threads" maxThreads="$threads" maxMemory=31 "$canu_input" "$r"; "
            assembly_commands+="date >> "$set_name"_canu_v1.8.time; "
            assembly_commands+="cp "$out_dir"/canu.contigs.fasta "$set_name"_canu_v1.8.fasta; "
            assembly_commands+="cp "$out_dir"/canu.contigs.gfa "$set_name"_canu_v1.8.gfa; "
            assembly_commands+="gzip "$set_name"_canu_v1.8.fasta "$set_name"_canu_v1.8.gfa; "
            assembly_commands+="rm -r "$out_dir";"
            sbatch --job-name=real_"$set_name"_canu --wrap "$assembly_commands" --time=4-00:00:00 --partition=comp --account md52 --qos=normal --ntasks=1 --mem=32000 --cpus-per-task="$threads" --output "$set_name"_canu_v1.8.out --error "$set_name"_canu_v1.8.err
            sleep 0.2

        done
    done
done
