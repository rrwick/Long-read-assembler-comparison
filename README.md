
# Benchmarking of long-read assemblers for prokaryote whole genome sequencing

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2702442.svg)](https://doi.org/10.5281/zenodo.2702442)


This repo contains the supplementary figures, scripts and data used for our paper comparing long-read assemblers:<br>
[Wick RR, Holt KE. Benchmarking of long-read assemblers for prokaryote whole genome sequencing. F1000Research. 2019;8(2138).](https://f1000research.com/articles/8-2138)

Are you interested in the older version of this comparison which was hosted here on GitHub? You can still find it [here](https://github.com/rrwick/Long-read-assembler-comparison/tree/96dfbe7e6edd6195cdd7fbe6f532f0022a7ebbb9).

<br><br>





## Figures

<p align="center"><img src="figures/Fig_1_simulated_read_sets.png" alt="Figure 1" width="90%"></p>

__Figure 1.__ Assembly results for the simulated read sets which cover a wide variety of parameters for length, depth and quality. 'Miniasm+' here refers to the entire Miniasm/Minipolish assembly pipeline. __A:__ Proportion of each possible assembly outcome. __B:__ Relative contiguity of the chromosome for each assembly, showing cleanliness of circularisation. __C:__ Relative contiguity of all plasmids in the assemblies, showing cleanliness of circularisation. __D:__ Sequence identity of each assembly's longest alignment to the chromosome. __E:__ The maximum indel error size in each assembly's longest alignment to the chromosome. __F:__ Total time taken (wall time) for each assembly. __G:__ Maximum RAM usage for each assembly.

<br><br><br><br>

<p align="center"><img src="figures/Fig_2_real_read_sets.png" alt="Figure 2" width="90%"></p>

__Figure 2.__ Assembly results for the real read sets, half containing ONT MinION reads (circles) and half PacBio RSII reads (X shapes). 'Miniasm+' here refers to the entire Miniasm/Minipolish assembly pipeline. __A:__ Proportion of each possible assembly outcome. __B:__ Relative contiguity of the chromosome for each assembly, showing cleanliness of circularisation. __C:__ Relative contiguity of all plasmids in the assemblies, showing cleanliness of circularisation. __D:__ Sequence identity of each assembly's longest alignment to the chromosome. __E:__ The maximum indel error size in each assembly's longest alignment to the chromosome. __F:__ Total time taken (wall time) for each assembly. __G:__ Maximum RAM usage for each assembly.

<br><br><br><br>





## Supplementary figures

<p align="center"><img src="figures/Fig_S01_ref_genome_replicons.png" alt="Figure S1" width="60%"></p>

__Figure S1.__ Distributions of chromosome sizes (__A__), plasmid sizes (__B__) and per-genome plasmid counts (__C__) for the reference genomes used to make the simulated read sets.

<br><br><br><br>



<p align="center"><img src="figures/Fig_S02_badread_parameters.png" alt="Figure S2" width="80%"></p>

__Figure S2.__ Badread parameter histograms for the simulated read sets. __A:__ Mean read depths were sampled from a uniform distribution ranging from 5x to 200x. __B:__ Mean read lengths were sampled from a uniform distribution ranging from 100 to 20000 bp. __C:__ Read length standard deviations were sampled from a uniform distribution ranging from 100 to twice that set's mean length (up to 40000 bp). __D:__ Mean read identities were sampled from a uniform distribution ranging from 80% to 99%. __E:__ Max read identities were sampled from a uniform distribution ranging from that set's mean identity plus 1% to 100%. __F:__ Read identity standard deviations were sampled from a uniform distribution ranging from 1% to the max identity minus the mean identity. __G, H and I:__ Junk, random and chimera rates were all sampled from an exponential distribution with a mean of 2%. __J:__ Glitch sizes/skips were sampled from a uniform distribution ranging from 0 to 100. __K:__ Glitch rates for each set were calculated from the size/skip according to this formula: 100000/(1.6986^(s/10)). __L:__ Adapter lengths were sampled from an exponential distribution with a mean of 50.

<br><br><br><br>



<p align="center"><img src="figures/Fig_S03_replicon_depths.png" alt="Figure S3" width="60%"></p>

__Figure S3.__ Top: the target simulated depth of each replicon relative to the chromosome. The smaller the plasmid, the wider the range of possible depths. Bottom: the absolute read depth of each replicon after read simulation.

<br><br><br><br>



<p align="center"><img src="figures/Fig_S04_commands.png" alt="Figure S4" width="80%"></p>

__Figure S4.__ Commands used for each of the eight assemblers tested.

<br><br><br><br>



<p align="center"><img src="figures/Fig_S05_circularisation.png" alt="Figure S5" width="80%"></p>

__Figure S5.__ Possible states for the assembly of a circular replicon. Reference sequences are shown in the inner circles in black and aligned contig sequences are shown in the outer circles in colour (red at the contig start to violet at the contig end). __A:__ Complete assembly with perfect circularisation. __B:__ Complete assembly but with missing bases leading to a gapped circularisation. __C:__ Complete assembly but with duplicated bases leading to overlapping circularisation. __D:__ Incomplete assembly due to fragmentation (multiple contigs per replicon). __E:__ Incomplete assembly due to missing sequence. __F:__ Incomplete assembly due to misassembly (non-contiguous sequence in the contig).

<br><br><br><br>



<p align="center"><img src="figures/Fig_S06_tripled_reference.png" alt="Figure S6" width="80%"></p>

__Figure S6.__ Reference triplication for assembly assessment. __A:__ Due to the ambiguous starting position of a circular replicon, a completely-assembled contig will typically not align to the reference in a single unbroken alignment. __B:__ Doubling the reference sequence will allow for a single alignment, regardless of starting position. __C:__ However, if the contig contains start/end overlap (i.e.\ contiguity >100%) then even a doubled reference may not be sufficient to achieve a single alignment, depending on the starting position. __D:__ A tripled reference allows for an unbroken alignment, regardless of starting position, even in cases of >100% contiguity.

<br><br><br><br>



<p align="center"><img src="figures/Fig_S07_problem_plots.png" alt="Figure S7" width="100%"></p>

__Figure S7.__ Contiguity of the simulated read set assemblies plotted against Badread parameters for each of the tested assemblers. These plots show how well the assemblers tolerate different problems in the read sets. __A:__ Mean read depth (higher is better). __B:__ Max read identity (higher is better). __C:__ N50 read length (higher is better). __D:__ The sum of random read rate and junk read rate (lower is better). __E:__ Chimeric read rate (lower is better). __F:__ Adapter sequence length (lower is better). __G:__ Glitch size/skip (lower is better).

<br><br><br><br>



<p align="center"><img src="figures/Fig_S08_simulated_plasmids.png" alt="Figure S8" width="75%"></p>

__Figure S8.__ Plasmid completion for the simulated read set assemblies for each of the tested assemblers, plotted with plasmid length and read depth. Solid dots indicate completely assembled plasmids (contiguity ≥99%) while open dots indicate incomplete plasmids (contiguity <99%). Percentages in the plot titles give the proportion of plasmids which were completely assembled.

<br><br><br><br>



<p align="center"><img src="figures/Fig_S09_real_plasmids.png" alt="Figure S9" width="75%"></p>

__Figure S9.__ Plasmid completion for the real read set assemblies for each of the tested assemblers, plotted with plasmid length and read depth. Solid dots indicate completely assembled plasmids (contiguity ≥99%) while open dots indicate incomplete plasmids (contiguity <99%). Percentages in the plot titles give the proportion of plasmids which were completely assembled.

<br><br><br><br>




## License

[Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/legalcode)

