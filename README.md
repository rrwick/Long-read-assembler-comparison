
# Benchmarking of long-read assemblers for prokaryote whole genome sequencing

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2702442.svg)](https://doi.org/10.5281/zenodo.2702442)


This repo contains the supplementary figures, scripts and data used for our upcoming paper comparing long-read assemblers. It should be available shortly – stay tuned!

Are you interested in the older version of this comparison which was hosted here on GitHub? You can still find it [here](https://github.com/rrwick/Long-read-assembler-comparison/tree/96dfbe7e6edd6195cdd7fbe6f532f0022a7ebbb9).

<br><br>



## Supplementary figures

<p align="center"><img src="supplementary_figures/Fig_S01_ref_genome_replicons.png" alt="Minipolish" width="60%"></p>

__Figure S1.__ Distributions of chromosome sizes (A), plasmid sizes (B) and per-genome plasmid counts (C) for the reference genomes used to make the simulated read sets.

<br><br><br><br>



<p align="center"><img src="supplementary_figures/Fig_S02_badread_parameters.png" alt="Minipolish" width="80%"></p>

__Figure S2.__ Badread parameter histograms for the simulated read sets. (A) Mean read depths were sampled from a uniform distribution ranging from 5x to 200x. (B) Mean read lengths were sampled from a uniform distribution ranging from 100 to 20000 bp. (C) Read length standard deviations were sampled from a uniform distribution ranging from 100 to twice that set's mean length (up to 40000 bp). (D) Mean read identities were sampled from a uniform distribution ranging from 80% to 99%. (E) Max read identities were sampled from a uniform distribution ranging from that set's mean identity plus 1% to 100%. (F) Read identity standard deviations were sampled from a uniform distribution ranging from 1% to the max identity minus the mean identity. (G, H and I) Junk, random and chimera rates were all sampled from an exponential distribution with a mean of 2%. (J) Glitch sizes/skips were sampled from a uniform distribution ranging from 0 to 100. (K) Glitch rates for each set were calculated from the size/skip according to this formula: 100000/(1.6986^(s/10)). (L) Adapter lengths were sampled from an exponential distribution with a mean of 50.

<br><br><br><br>



<p align="center"><img src="supplementary_figures/Fig_S03_replicon_depths.png" alt="Minipolish" width="60%"></p>

__Figure S3.__ Top: the target simulated depth of each replicon relative to the chromosome. The smaller the plasmid, the wider the range of possible depths. Bottom: the absolute read set of each replicon after read simulation.

<br><br><br><br>



<p align="center"><img src="supplementary_figures/Fig_S04_commands.png" alt="Minipolish" width="80%"></p>

__Figure S4.__ Commands used for each of the six assemblers tested.

<br><br><br><br>



<p align="center"><img src="supplementary_figures/Fig_S05_circularisation.png" alt="Minipolish" width="80%"></p>

__Figure S5.__ Possible states for the assembly of a circular replicon. Reference sequences are shown in the inner circles in black and aligned contig sequences are shown in the outer circles in colour (red at the contig start to violet at the contig end). (A) Complete assembly with perfect circularisation. (B) Complete assembly but with missing bases leading to a gapped circularisation. (C) Complete assembly but with duplicated bases leading to overlapping circularisation. (D) Incomplete assembly due to fragmentation (multiple contigs per replicon). (E) Incomplete assembly due to missing sequence. (F) Incomplete assembly due to misassembly (non-contiguous sequence in the contig).

<br><br><br><br>



<p align="center"><img src="supplementary_figures/Fig_S06_tripled_reference.png" alt="Minipolish" width="80%"></p>

__Figure S6.__ Reference triplication for assembly assessment. (A) Due to the ambiguous starting position of a circular replicon, a completely-assembled contig will typically not align to the reference in a single unbroken alignment. (B) Doubling the reference sequence will allow for a single alignment, regardless of starting position. (C) However, if the contig contains start/end overlap (i.e.\ contiguity >100%) then even a doubled reference may not be sufficient to achieve a single alignment, depending on the starting position. (D) A tripled reference allows for an unbroken alignment, regardless of starting position, even in cases of >100% contiguity.

<br><br><br><br>



<p align="center"><img src="supplementary_figures/Fig_S07_problem_plots.png" alt="Minipolish" width="100%"></p>

__Figure S7.__ Contiguity of the simulated read set assemblies plotted against Badread parameters for each of the tested assemblers. These plots show how well the assemblers tolerate different problems in the read sets. (A) Mean read depth (higher is better). (B) Max read identity (higher is better). (C) N50 read length (higher is better). (D) The sum of random read rate and junk read rate (lower is better). (E) Chimeric read rate (lower is better). (F) Adapter sequence length (lower is better). (G) Glitch size/skip (lower is better).

<br><br><br><br>



<p align="center"><img src="supplementary_figures/Fig_S08_simulated_plasmids.png" alt="Minipolish" width="75%"></p>

__Figure S8.__ Plasmid completion for the simulated read set assemblies for each of the tested assemblers, plotted with plasmid length and read depth. Solid dots indicate completely assembled plasmids (contiguity ≥99%) while open dots indicate incomplete plasmids (contiguity <99%). Percentages in the plot titles give the proportion of plasmids which were completely assembled.

<br><br><br><br>



<p align="center"><img src="supplementary_figures/Fig_S09_real_plasmids.png" alt="Minipolish" width="75%"></p>

__Figure S9.__ Plasmid completion for the real read set assemblies for each of the tested assemblers, plotted with plasmid length and read depth. Solid dots indicate completely assembled plasmids (contiguity ≥99%) while open dots indicate incomplete plasmids (contiguity <99%). Percentages in the plot titles give the proportion of plasmids which were completely assembled.

<br><br><br><br> 



<p align="center"><img src="supplementary_figures/Fig_S10_plasmid_contiguity.png" alt="Minipolish" width="100%"></p>

__Figure S10.__ The relative contiguity of the plasmids for each real read set assembly (A) and simulated read set assembly (B).

<br><br><br><br>




## License

[Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/legalcode)

