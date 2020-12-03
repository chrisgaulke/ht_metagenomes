# ht_metagenomes
Analysis of the performance of different metagenomic library prep protocols
 when applied to variable community types and input concentrations. 

---
---

## Structure 

There are two main subdirectories: 1) analysis and 2) data. The directory 
analysis/ contains all of the combined data and metadata files needed for the 
analysis, while data/ houses the original data files that were used to generate
these files

### Code files

These files are found in the top level directory and can be executed in any 
order with the exception of combine_data.R, which must be run first to generate 
the input files. The files are: 

1. combine_data.R
- combines individual input data files into one data frame and writes out these 
data frames
2. ht_metagenomes_seq_stats.R
-conducts stats on library and seqeunce quality
3. 2020_03_25_ht_metagenomes_genes.R
-conducts gene level analyses
4. ht_meta_tax_revised_2020_02_24.R
-conducts metaphlan vs mock analysis
5. 2020_01_16_kraken2_mock_analysis.R
-conducts kraken2 vs mock analysis
6. 2020_03_23_ht_metagenomes_metaphlan.R
- conducts metaphlan taxonomic analysis 
7. 2020_04_02_ht_metagenomes_kraken2.R
- conducts kraken taxonomic analysis 


