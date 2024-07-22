# Pogoniulus-phylogenomics-manuscript-data

This repository contains supplementary information for the manuscript "Introgression across narrow contact zones shapes the genomic landscape of phylogenetic variation in a tinkerbird clade", as well as scripts and data necessary to reproduce analyses and figures.
Files are organized as follow:
1. Appendix1_ddRADseq_samples_metadata.xlsx: A sreapsheet giving the geographical origin, distance to contact zone and Red-fronted tinkerbird ancestry (estimated with Admixture) for each of the samples included in the ddRADseq datasets. ddRADseq_samples_metadata_South_Africa.csv and ddRADseq_samples_metadata_Uganda-Kenya.csv include the same information in csv format for easy importation in R.
2. Appendix2_coordinates_mt_genes.xlsx: A spreasheet listing the 150 mito-nuclear genes investigated, with their coordinates in the Pogoniulus pusillus reference genome.
3. Linkage_map_YFT_LDhat_100kb.txt: The linkage map of Pogoniulus extoni, inferred with LDhat.
4. Neighbour_Joining_trees: A directory including Neighbour-Joining trees calculated in sliding windows using the Whole-Genome Sequencing data. Trees were calulated for each chromosome separately (SUPER_1 to SUPER_44 and SUPER_Z). For each chromosome the following files are included:
   * SUPER_*_all_500S_maxmiss60_NJ_trees.trees: a files containing the trees calculated for this chromosome on separate lines, or NA if NJ calculation failed.
   * SUPER_*_all_500S_maxmiss60_windows_stats.tsv: a table with metadata for the genomic windows from which the NJ trees were calculated. Each line in this table corresponds to the same line number in the previous file.
   * SUPER_*_all_500S_maxmiss60_NJ_ASTRAL.trees: the trees, filtered to use as input for ASTRAL.
   * SUPER_*_all_500S_maxmiss60_NJ_TWISST.trees: the trees, filtered to use as input for TWISST.
     
Each of these files is provided for four sample subsets: all samples (SUPER_*_all), core range samples (SUPER_*_allopatric), contact zones samples (SUPER_*_sympatric) and control samples (SUPER_*_control).

5. TWISST_weights: A directory including topology weights for each chromosome, as output by TWISST. Again, four files are provided for each chromosome, corresponding to the four samples subsets.
6. HZAR_single_SNPs_South_Africa and HZAR_single_SNPs_Uganda-Kenya: Two directories containing the results of single SNP HZAR cline analyses in the two contact zones. The analysis was run on batches of 50 SNPs, and the results are given in the files "Estimates_hzar_*.txt, which are tables where each line is a SNP and the columns are the estimated cline parameters. These directories also contain the HZAR input files (Allele_frequencies*HZAR.csv) generated with Stacks populations, and the coordinates of the SNPs in plink format (SNPs*.plink.map).
7. Mito-nuclear_genes: A directory containing files necessary to repeat the analyses of mito-nuclear genes.

In addition to that, the R scripts used to analyse the data and produce figures are given. They can be executed from the root of the directory. (Note: the scripts were written and run in R v. 4.1.1). 
