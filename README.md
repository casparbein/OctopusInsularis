# Oceanic currents and habitat breaks as drivers of the diversification of Octopus insularis

In this repository, I lay out the worflow I followed to analyse ddRADseq data for the project "Oceanic currents and habitat breaks as drivers of the diversification of Octopus insularis" as specified in the Materials/Methods section of the manuscript. For each part of the analyes, I provide relevant scripts in this repository.

## Basic workflow

Part of the analyses has been conducted on a remote machine via SLURM scripts. As these scripts are not universally applicable, I provide an example script to demonstrate its overall appearance and present the commands I used to run the analysis here.

- relevant scripts in this repository: [example.sh](/example.sh/)


### ddRADseq assembly with iPyrad

I demulitplexed and assembled raw reads of 71 individuals that had been sequenced using the double-digest restriction association DNA sequencing (ddRADseq, 100 bp single end reads) approach using ipyrad v.0.9.50 . For more information on how to proceed with ipyrad, please check https://ipyrad.readthedocs.io/en/master/. I reassembled the data after removing two individuals with low read numbers for a final assembly of 69 individuals.

In general, I followed the usual seven-steps approach: `ipyrad -p params-_90_69.txt -s 234567 -t 1 -c 20 -f`.

As run times for this assembly routinely exceeded 48 hours and I had to override unfinished assemblies, I used the `-f` and `-c 20` flags to assign 20 cores and force overrides. In the parameter file, I kept the default values for all but the sequence similarity threshold, `[clust_threshold]` (I tried 0.9 and 0.95) and the minimum number of individuals a locus had to be found in to be considered `[min_samples_locus]` to 4.

- relevant scripts in this repository: [ipyrad](/ipyrad/)

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

- relevant scripts in this repository: [filtering](/filtering/)

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

- relevant scripts in this repository: [EMU PCA](/EMU_PCA/)

#### tess3R

Another approach to infer population structure I used is tess3R. Similar to well-known clustering algorithms like STRUCTURE or ADMIXTURE, this approach tries assign individuals to K ancestral clusters. 
Here, genotype matrix (at each biallelic locus, a SNP can be translated into 0: homozygous for allelle 1, 1: heterozygous, 2: homozygous for allele 2) is factorized into a matrix Q of size n (number of individuals) x K (number of clusters), describing the assignment probability of any given individual i to cluster Kj, and another matrix G with size K x g (number of loci), describing the ancestral allele frequencies of gx in cluster Kj. Additionally, it uses the geographic location of sampled individuals as a prior for the assignment of individuals to clusters, and to perform a spatial interpolation that reflects the most likely assignment in unsampled areas.
I used K=1 to K=11 (the number of populations), similar to the eigenvalues in the EMU-PCA to cover a broad range of possible clusters and used both the more permissive ("69inds_40MD") and more stringent ("64inds_20MD") datasets for inference.

For more information on tess3R, see: https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12471 . The program is implemented as an R package and the documentation can be found here: https://bcm-uga.github.io/TESS3_encho_sen/index.html .

- relevant scripts in this repository: [tess3R](/tess3R/)

#### fineRADstructure

Another clustering algorithm specifically designed for RADseq data is fineRADstructure. By calculating co-ancestry between all pairs of individuals (using all SNPs within one locus), it can detect hierarchical population structure without a prior based on sampling location. Since fineRADstructure can also handle large amounts of missing data, I ran the program with the original dataset as well as with the "69inds_40MD" and "64inds_20MD" datasets. First, I converted the vcf files to the domestic fineRADstructure format: 
`/opt/miniconda3/bin/RADpainter hapsFromVCF name.vcf`. 

Then, I inferred clusters of individuals based on shared co-ancestry: 
`/opt/miniconda3/bin/finestructure -x 100000 -y 1000000 -z 1000 name_chunks.out  name_mcmc.xml`. 

Lastly, I computed a clustering tree: 
`opt/miniconda3/bin/finestructure -m T -x 100000 name_chunks.out  name_mcmc.xml name_TREEmcmc.xml`. 

A manual of how to use fineRADstructure can be found here: https://www.milan-malinsky.org/fineradstructure and R scripts for plotting the results were retrieved from the developer's github page: https://github.com/millanek/fineRADstructure. The original publication can be read here: https://academic.oup.com/mbe/article/35/5/1284/4883220?login=false .

- relevant scripts in this repository: [fineRADstructure](/fineRADstructure/)

### Phylogenetic analysis

#### IQtree

To investigate maximum likelihood phylogenetic relationships, I used IQtree v.2.1.2. The program is straightforward to use and includes options for model detection (ModelFinder), ultrafast bootstraps and replicated runs. I ran ML inference for both the more permissive ("69inds_40MD") and the ("64inds_20MD") datasets, which I converted to the fasta file format with PGDSpider, five times independently (`--runs 5`) with 1000 ultrafast bootstraps (`-B 1000`) and a maximum of 1500 iterations (`-nm 1500`):

`/opt/bin/iqtree2 -s input.fa -m MFP -B 1000 -T AUTO -ntmax 8 --seqtype DNA -nm 1500 --runs 5 --seed 35681 `

A comprehensive manual including all relevant publications of IQtree is available here: http://www.iqtree.org/.
I plotted selected trees both with an arbitrary root (outgroup S-Coastal and S-Oceanic) and unrooted with ggtree in R. The extensive ggtree manual can be found here:https://yulab-smu.top/treedata-book/index.html, and the original publication here:https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12628 .
.

- relevant scripts in this repository: [IQtree](/IQtree/)

#### SVDquartets / tetrad

A coaslescent based approach I used for tree reconstruction is SVDquartets (Singular Value Decomposition quartets) as implemented in the tetrad application in ipyrad. The SVDquartets algorithm computes true splits of any combination of four individuals in the dataset into two sets based on unlinked SNPs that are summarized in a flattened tensor of size 16 x 16. 
Tetrad randomly chooses one SNP per locus for tree inference and bootstrapping. This approach maximizes the amount of SNPs used in the analysis, and bootstrapped trees hold additional information not visible in the primarily inferred tree. 
I ran the SVDquartets analysis in tetrad for the "full", "69inds_40MD" and "64inds_20MD" datasets within a jupyter notebook, using 100 non-parametric bootstrap replicates for each dataset. 

The tetrad manual is available here: https://ipyrad.readthedocs.io/en/master/API-analysis/cookbook-tetrad.html, and the original publication can be found here: https://academic.oup.com/bioinformatics/article/30/23/3317/206559?login=false .

- relevant scripts in this repository: [tetrad](/tetrad_analysis.ipynb)

### Population summary statistics

I calculated Fit (individual inbreeding coefficient) and individual missingness for the "69inds_40MD" and "64inds_20MD" datasets with vcftools v.0.1.14 like this:
- individual inbreeding coefficient:
    `vcftools --vcf name.vcf --het --out name_full_het`     
- individual missingness:
    `vcftools --vcf name.vcf --missing-indv --out name_miss`   

In R, I tested if differences in Fit and pi SNP (calculated with DnaSP) of individuals belonging to different clusters (as defined by the population structure analyses) were significant with a BH-corrected pairwise t-test and created plots of Fit and pi SNP stratified by cluster assigment. I also performed a linear regression to test if missingness was significantly correlated to Fit and piSNP (result: p>0.05 in both cases).

- relevant scripts in this repository: [population summary statistics](/pop_summary_statistics/)

### Demographic history

I estimated past changes in effective population size for N-Coastal, S-Coastal and S-Oceanic and isolation/geneflow between N-Coastal/S-Coastal and S-Oceanic/S-Coastal with the program dadi. It uses a diffusion approximation approach to fit parameters of predefined models to the observed site frequency spectrum (SFS) by simulating a SFS and optimizing fit between simulated and observed SFS. Dadi does not work well with linked SNPs and large amounts of missing data, so perfiltering of the "full" vcf file was obligatory.
To maximize to number of shared SNPs within a certain population by simultaneously minimizing SNP linkage, I created costumized vcf input files for each analysis using vcftools.

Single population modelling:

- N-Coastal: `vcftools --vcf ../ipyrad_assembly/clean_90_69_real_outfiles/name.vcf --max-missing 0.6 --thin 1000 --remove remove_all_but_CEN --out CEN_str_full_06_thin --recode`

- S-Coastal: `vcftools --vcf ../ipyrad_assembly/clean_90_69_real_outfiles/name.vcf --max-missing 0.6 --thin 1000 --remove remove_all_but_ALBA --out ALBA_str_full_06_thin --recode`

- S-Oceanic: `vcftools --vcf ../ipyrad_assembly/clean_90_69_real_outfiles/name.vcf --max-missing 0.6 --thin 1000 --remove remove_all_but_TM --out TM_str_full_06_thin --recode`

Two population modelling: 
- N-Coastal/S-Coastal: `vcftools --vcf ../ipyrad_assembly/clean_90_69_real_outfiles/name.vcf --max-missing 0.6 --thin 1000 --remove remove_all_but_ALBA_CEN --out ALBACEN_str_full_06_thin --recode`

- S-Oceanic/S-Coastal: `vcftools --vcf ../ipyrad_assembly/clean_90_69_real_outfiles/name.vcf --max-missing 0.6 --thin 1000 --remove remove_all_but_ALBA_TM --out TMALBA_str_full_06_thin --recode`

To retain as many loci as possible while removing as much missing data as possible, I down-projected all input vcf files. For this purpose, as well as for streamlining and automizing the model fitting process, I used Daniel Portik`s dadi pipeline. All models I included in the analysis are implemented in the dadi pipeline, too. For the two population modelling, I performed a goodness-of-fit test (also part of the dadi pipeline), where the program simulates an SFS based on the optimized parameters obtained from model fitting, and runs optimizations on this simulated SFS. If the previously fit model is a "good fit" to the data, the original log. likelihood value will be included in the distribution of log. likelihood values from fits to the simulated SFS.

For more information on how to use dadi, please refer to the manual here: https://dadi.readthedocs.io/en/latest/ , and the original publication here: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000695 . The dadi pipeline can be found on github:https: https://github.com/dportik/dadi_pipeline , and the original publication is available here: https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.14266 .

- relevant scripts in this repository: [dadi](/dadi/)

### References

- Avni, E., Cohen, R., & Snir, S. (2015). Weighted quartets phylogenetics. Systematic biology, 64(2), 233-242.

- Barratt, C. D., Bwong, B. A., Jehle, R., Liedtke, H. C., Nagel, P., Onstein, R. E., ... & Loader, S. P. (2018). Vanishing refuge? Testing the forest refuge hypothesis in coastal East Africa using genome‐wide sequence data for seven amphibians. Molecular Ecology, 27(21), 4289-4308.

- Bates, D., Sarkar, D., Bates, M. D., & Matrix, L. (2007). The lme4 package. R package version, 2(1), 74.

- Caye, K., Deist, T. M., Martins, H., Michel, O., & François, O. (2016). TESS3: fast inference of spatial population structure and genome scans for selection. Molecular Ecology Resources, 16(2), 540-548.

- Chifman, J., & Kubatko, L. (2014). Quartet inference from SNP data under the coalescent model. Bioinformatics, 30(23), 3317-3324.

- Danecek, P., Auton, A., Abecasis, G., Albers, C. A., Banks, E., DePristo, M. A., ... & 1000 Genomes Project Analysis Group. (2011). The variant call format and VCFtools. Bioinformatics, 27(15), 2156-2158.

- Eaton, D. A., & Overcast, I. (2020). ipyrad: Interactive assembly and analysis of RADseq datasets. Bioinformatics, 36(8), 2592-2594.

- Gutenkunst, R. N., Hernandez, R. D., Williamson, S. H., & Bustamante, C. D. (2009). Inferring the joint demographic history of multiple populations from multidimensional SNP frequency data. PLoS genet, 5(10), e1000695.

- Hoang, D. T., Chernomor, O., Von Haeseler, A., Minh, B. Q., & Vinh, L. S. (2018). UFBoot2: improving the ultrafast bootstrap approximation. Molecular biology and evolution, 35(2), 518-522.

- Kalyaanamoorthy, S., Minh, B. Q., Wong, T. K., Von Haeseler, A., & Jermiin, L. S. (2017). ModelFinder: fast model selection for accurate phylogenetic estimates. Nature methods, 14(6), 587-589.

- Malinsky, M., Trucchi, E., Lawson, D. J., & Falush, D. (2018). RADpainter and fineRADstructure: population inference from RADseq data. Molecular biology and evolution, 35(5), 1284-1290.

- Meisner, J., Liu, S., Huang, M., & Albrechtsen, A. (2021). Large-scale inference of population structure in presence of missingness using PCA. Bioinformatics, 37(13), 1868-1875.

- Minh, B. Q., Schmidt, H. A., Chernomor, O., Schrempf, D., Woodhams, M. D., Von Haeseler, A., & Lanfear, R. (2020). IQ-TREE 2: new models and efficient methods for phylogenetic inference in the genomic era. Molecular biology and evolution, 37(5), 1530-1534.

- Portik, D. M., Leaché, A. D., Rivera, D., Barej, M. F., Burger, M., Hirschfeld, M., ... & Fujita, M. K. (2017). Evaluating mechanisms of diversification in a Guineo‐Congolian tropical forest frog using demographic model selection. Molecular ecology, 26(19), 5245-5263.

- Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M. A., Bender, D., ... & Sham, P. C. (2007). PLINK: a tool set for whole-genome association and population-based linkage analyses. The American journal of human genetics, 81(3), 559-575.

- Yu, G., Smith, D. K., Zhu, H., Guan, Y., & Lam, T. T. Y. (2017). ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods in Ecology and Evolution, 8(1), 28-36.
