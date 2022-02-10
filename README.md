# Oceanic currents and habitat breaks as drivers of the diversification of Octopus insularis

In this repository, I have uploaded all scripts I used to analyse ddRADseq data for the project "Oceanic currents and habitat breaks as drivers of the diversification of Octopus insularis" as specified in the Materials/Methods section of the manuscript. For each part of the analyes, relevant scripts can be found in the respective folders of this repository.

## Basic workflow

Part of the analyses has been conducted on a remote machine via SLURM scripts. As these scripts are not universally applicable, I provide an example script to demonstrate its overall appearance and present the commands I used to run the analysis here.


### ddRADseq assembly with iPyrad

I demulitplexed and assembled raw reads of 71 individuals that had been sequenced using the double-digest restriction association DNA sequencing (ddRADseq, 100 bp single end reads) approach using ipyrad v.0.9.50 . For more information on how to proceed with ipyrad, please check https://ipyrad.readthedocs.io/en/master/. I reassembled the data after removing two individuals with low read numbers for a final assembly of 69 individuals.

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

Therefore, I decided to use another PCA approach called EMU-PCA (Expectation Maximization-PCA for Ultra-low Coverage Sequencing Data), which was designed to perform missing value imputations not by drawing from a predefined allele frequency distribution, but by Singular Value Decomposition (SVD) and low-rank approximation. Details about the algorithm can be found here: https://academic.oup.com/bioinformatics/article-abstract/37/13/1868/6103565?redirectedFrom=fulltext, and the program is available here: https://github.com/Rosemeis/emu.

EMU PCA uses .bed files (three files with extensions .bed, .fam, .bed) as input, so I had to convert my vcf input files "69inds_40MD" and "64inds_20MD" data to this format with plink v. v1.90b6.22. This can be tricky with .vcf files produced by iyprad, since every RAD locus has its own #CHROM name. Plink on the other hand only accepts around 10,000 unique chromosome names and throws an error otherwise. To circumvent this, I renamed each #CHROM to 0 (any other name should work, too) and reassigned the position of each SNP by adding the position of the previous SNP to its position. This is necessary because many positions will be the same across loci and will unevitably result in conflicting double-assignments when present on a single artifical chromosome. 

I then converted the manipulated .vcf file to .bed format using the following command:

`plink   --allow-extra-chr --make-bed --out 69inds_40MD --recode --vcf 69inds_40MD.recode.vcf`

The `--allow extra-chr` flag was not necessary in this case, but will be necessary if there are more than 1 unique #CHROM name in the .vcf file.

With EMU-PCA, to retrieve PC scores for all individuals depending on different numbers of eigenvalues, I employed the following command:

`python ../emu/emu.py -p 69inds_40MD -t 10 --iter 500 -o 69inds_40MD_2 -f 0.001 -e X`

Since eigenvalues in the EMU-PCA can be interpreted as K clusters in any structuring algorithm, I tested eigenvalues from X=2 to X=11, which in turn will also embody the rank of the matrix that is used to impute missing values. The analysis failed if I did not filter out singletons, so I included the `-f` flag.

I include an example output file of an EMU-PCA run, as well as the python script with which I prepared .vcf files for plink, in this repository. 

#### tess3R

Another approach to infer population structure I used is tess3R. Similar to well-known clustering algorithms like STRUCTURE or ADMIXTURE, this approach tries assign individuals to K ancestral clusters. 
Here, genotype matrix (at each biallelic locus, a SNP can be translated into 0: homozygous for allelle 1, 1: heterozygous, 2: homozygous for allele 2) is factorized into a matrix Q of size n (number of individuals) x K (number of clusters), describing the assignment probability of any given individual i to cluster Kj, and another matrix G with size K x g (number of loci), describing the ancestral allele frequencies of gx in cluster Kj. Additionally, it uses the geographic location of sampled individuals as a prior for the assignment of individuals to clusters, and to perform a spatial interpolation that reflects the most likely assignment in unsampled areas.
I used K=1 to K=11 (the number of populations), similar to the eigenvalues in the EMU-PCA to cover a broad range of possible clusters and used both the more permissive ("69inds_40MD") and more stringent ("64inds_20MD") datasets for inference.

For more information on tess3R, see: https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12471 . The program is implemented as an R package and the documentation can be found here: https://bcm-uga.github.io/TESS3_encho_sen/index.html .

The R script that I used to perform the tess3R analysis is included in this repository.

#### fineRADstructure
TBC.





 

