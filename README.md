# Oceanic currents and habitat breaks as drivers of the diversification of Octopus insularis

In this repository, I have uploaded all scripts I used to analyse ddRADseq data for the project "Oceanic currents and habitat breaks as drivers of the diversification of Octopus insularis" as specified in the Materials/Methods section of the manuscript. For each part of the analyes, relevant scripts can be found in the respective folders of this repository.

## Basic workflow

Part of the analyses has been conducted on a remote machine via SLURM scripts. As these scripts are not universally applicable, I provide an example script to demonstrate its overall appearance and present the commands I used to run the analysis here.


### ddRADseq assembly with iPyrad

I demulitplexed and assembled raw reads of 71 individuals that had been sequenced using the double-digest restriction association DNA sequencing (ddRADseq) approach using ipyrad v.0.9.50. For more information on how to proceed with ipyrad, please check https://ipyrad.readthedocs.io/en/master/. I reassembled the data after removing two individuals with low read numbers for a final assembly of 69 individuals.

In general, I followed the usual seven-steps approach: `ipyrad -p params-_90_69.txt -s 234567 -t 1 -c 20 -f`
As run times for this assembly routinely exceeded 48 hours and I had to override unfinished assemblies, I used the `-f` and `-c 20` flags to assign 20 cores and force overrides. In the parameter file, I kept the default values for all but the sequence similarity threshold, `[clust_threshold]` (I tried 0.9 and 0.95) and the minimum number of individuals a locus had to be found in to be considered `[min_samples_locus]` to 4.

One can find the parameters file for the 0.9 sequence similarity threshold and all summary statistics for both 0.9 and 0.95 sequence similarity thresholds (in this case including the individuals with an insufficient amount of reads) in this repository.

### Filtering for missing data

After assembly, I used the .loci and .vcf output file formats for downstream analysis. Expressed as a single nucleotide polymorphism (SNP) matrix, the .vcf file contained 66.34 % missing data. Therefore, I filtered the raw vcf file using vcftools v.0.1.14 for downstream analysis creating several datasets as follows:
- a more permissive dataset filtered for a maximum of 40% missing data across the SNP matrix (69inds_40MD) 

    `vcftools --vcf ../ipyrad_assembly/clean_90_69_outfiles/clean_90_69.vcf --max-missing 0.6 --out 69inds_40MD --recode`
 
- a more stringent dataset filtered for a maximum of 20% missing data across the SNP matrix, also removing the five individuals with the highest amount of missing data (64inds_20MD)

    `vcftools --vcf ../ipyrad_assembly/clean_90_69_outfiles/clean_90_69.vcf --remove remove_inds --out 64inds_full --recode`
    
    `vcftools --vcf ../64inds_full/64inds_full.recode.vcf --max-missing 0.8 --out 64inds_20MD --recode`
    
    removed individuals were: BA18, OIC1, RN13, STH2, RN2A
    
 - various datasets filtered for 40% missing data with the 5 individuals with most missing data removed and thinned out to one SNP per locus for a subsequent demographic modelling analysis with dadi (see section about dadi for more details)

    `vcftools --vcf ../ipyrad_assembly/clean_90_69_outfiles/clean_90_69.vcf --max-missing 0.6 --thin 1000 --remove remove_necessary --out population_str_full_06_thin --recode`

A summary script of the vcftools filtering in plain text format can be found in this repository.

### Population structure analysis

#### EMU PCA

Even after filtering, our SNP matrix still contained a considerable amount of missing data. For any kind of principal component analysis (PCA), missing data has to be imputed one way or another. In iyprad, imputation of missing values for PCA can be done either based on pre-assigned groupings, where missing alleles at a certain locus are drawn from a distribution of allele frequencies on that locus in the respective grouping, or based on the whole sample, where alleles are then drawn from sample-wide frequencies. Both approaches can introduce biases on PC clustering, as samples with a lot of missing values tend to cluster with samples from preassigned groups or with each other, respectively.

Therefore, I decided to use another PCA approach called EMU-PCA (Expectation Maximization-PCA for Ultra-low Coverage Sequencing Data), which was designed to perform missing value imputations not by drawing from a predefined allele frequency distribution, but by Singular Value Decomposition (SVD) and low-rank approximation. Details about the algorithm can be found here: https://academic.oup.com/bioinformatics/article-abstract/37/13/1868/6103565?redirectedFrom=fulltext, and a handbook of how to use EMU PCA is available here: https://github.com/Rosemeis/emu.

TBC.


 

