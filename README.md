# stochastic exclusive gene choice
This set of software is to accompany the article (publication pending), please refer to it for more details on the project itself.

- [1. Requirements/Dependencies](#1-requirementsdependencies)
- [2. Introduction](#2-introduction)
- [3. Important notes](#3-important-notes)
  - [3.1. Nextflow caching and storage capacity](#31-nextflow-caching-and-storage-capacity)
  - [3.2. Results directory structure](#32-results-directory-structure)
- [4. Software description](#4-software-description)
  - [4.1. Data Preparation](#41-data-preparation)
    - [4.1.1. Prepare family set (PANTHER): **preprocess_panther.py**](#411-prepare-family-set-panther-preprocess_pantherpy)
    - [4.1.2. Prepare TPM datasets: **clean_dataset.ipynb**](#412-prepare-tpm-datasets-clean_datasetipynb)
  - [4.2. Chromosome-wise pipelines](#42-chromosome-wise-pipelines)
    - [4.2.1. Filter chromosomes: **bin/split_chromosomes.py**](#421-filter-chromosomes-binsplit_chromosomespy)
    - [4.2.2. Calculate thresholds and dichotomise: **bin/dichotomise.py**](#422-calculate-thresholds-and-dichotomise-bindichotomisepy)
    - [4.2.3. Filter genes: **bin/filter_genes.py**](#423-filter-genes-binfilter_genespy)
    - [4.2.4. Calculate IC: **bin/calculate_IC.py**](#424-calculate-ic-bincalculate_icpy)
    - [4.2.5. Calculate IC for shuffled chromosomes: **bin/calculate_chromosomes_shuffled_IC.py**](#425-calculate-ic-for-shuffled-chromosomes-bincalculate_chromosomes_shuffled_icpy)
    - [4.2.6. Collect shuffled IC: **bin/collect_chromosomes_shuffled_IC.py**](#426-collect-shuffled-ic-bincollect_chromosomes_shuffled_icpy)
  - [4.3. Family-wise pipeline](#43-family-wise-pipeline)
    - [4.3.1. Transform family set: **bin/transform_familyset.py**](#431-transform-family-set-bintransform_familysetpy)
    - [4.3.2. Dichotomise: **bin/dichotomise_families.py**](#432-dichotomise-bindichotomise_familiespy)
    - [4.3.3. Calculate IC: **bin/calculate_family_IC.py**](#433-calculate-ic-bincalculate_family_icpy)
    - [4.3.4. Permute families: **bin/permute_families.py**](#434-permute-families-binpermute_familiespy)
    - [4.3.5. Calculate shuffled IC: **bin/calculate_family_IC.py**](#435-calculate-shuffled-ic-bincalculate_family_icpy)
    - [4.3.6. Collect Shuffled IC: **bin/collect_shuffled_IC.py**](#436-collect-shuffled-ic-bincollect_shuffled_icpy)
  - [4.4. Chromosome-wise analysis](#44-chromosome-wise-analysis)
    - [4.4.1. Statistical tests: **permutation_test.ipynb**](#441-statistical-tests-permutation_testipynb)
    - [4.4.2. Plot statistical tests: **plot_chr_stat_tests.ipynb**](#442-plot-statistical-tests-plot_chr_stat_testsipynb)
    - [4.4.3. Plot thresholds/dichotomisation: **plot_thresholds.ipynb**](#443-plot-thresholdsdichotomisation-plot_thresholdsipynb)
    - [4.4.4. Overlap of tails of chromosomal stretches: **chromosome_overlap_analysis.ipynb**](#444-overlap-of-tails-of-chromosomal-stretches-chromosome_overlap_analysisipynb)
  - [4.5. Family-wise analysis](#45-family-wise-analysis)
    - [4.5.1. Statistical tests: **family_wise_analysis.ipynb**](#451-statistical-tests-family_wise_analysisipynb)
    - [4.5.2. Overlap of tails of families: **family_overlap_analysis.ipynb**](#452-overlap-of-tails-of-families-family_overlap_analysisipynb)
    - [4.5.3. Plot statistical tests: **plot_family_stat_tests.ipynb**](#453-plot-statistical-tests-plot_family_stat_testsipynb)
    - [4.5.4. Plot chosen families: **plot_family_overlap_heatmap.ipynb**](#454-plot-chosen-families-plot_family_overlap_heatmapipynb)
  - [4.6. Both chromosome- and family-wise analysis](#46-both-chromosome--and-family-wise-analysis)
    - [4.6.1. Cell type overlaps: **dataset_overlaps.ipynb**](#461-cell-type-overlaps-dataset_overlapsipynb)
  - [4.7. Misc](#47-misc)
    - [4.7.1. Transform Fitted Distributions data (from Mathematica): **fitted_distributions.ipynb**](#471-transform-fitted-distributions-data-from-mathematica-fitted_distributionsipynb)
    - [4.7.2. Dichotomisation Correlation: **dichotomisation_correlation.ipynb**](#472-dichotomisation-correlation-dichotomisation_correlationipynb)
    - [4.7.3. Thresholds for UMI (no bimodality): **alveolar_exploratory.ipynb**](#473-thresholds-for-umi-no-bimodality-alveolar_exploratoryipynb)
    - [4.7.4. Distribution fitting: **mathematica_code/Dichotomisation.nb**](#474-distribution-fitting-mathematica_codedichotomisationnb)
    - [4.7.5. Permutation test: **mathematics_code/permutation_test.nb**](#475-permutation-test-mathematics_codepermutation_testnb)
  - [4.8. Utility functions](#48-utility-functions)
    - [4.8.1. bin/pipeline_utils/gff_utils.py](#481-binpipeline_utilsgff_utilspy)
    - [4.8.2. bin/pipeline_utils/ic_utils.py](#482-binpipeline_utilsic_utilspy)
    - [4.8.3. bin/pipeline_utils/nextflow_utils.py](#483-binpipeline_utilsnextflow_utilspy)
    - [4.8.4. bin/pipeline_utils/pomegranate_parser.py](#484-binpipeline_utilspomegranate_parserpy)
    - [4.8.5. bin/pipeline_utils/vrs_utils.py](#485-binpipeline_utilsvrs_utilspy)

# 1. Requirements/Dependencies

* Standard Anaconda installation (pandas < 1.0, numpy, Jupyter notebooks, etc).
* tqdm -- for progress bars
* swifter -- for parallelisation (speeds the dichotomisation up drastically)
* pomegranate -- library for probability distributions
* Nextflow -- for pipelines [(link)](https://www.nextflow.io/)

The conda environment file (**stochastic_env.yml**) is provided, and it is encouraged to use it, so no library version mismatch happen.

# 2. Introduction

Different parts of the project were implemented as either standalone python scripts, scripts as the part of the nextflow pipelines, Mathematica notebooks, or Jupyter notebooks (mostly for the later stages of the project)

Note: even though the pipeline scripts can be run manually, I highly suggest running them via nextflow (this allows running the whole pipeline from input to the final results without interruption and any intervention from the user)

The following documentation describes all the components in the roughly logical order (more or less in the order parts of the project are to be executed).

For the details on the methods and the project itself see the article: [Gene Families With Stochastic Exclusive Gene Choice Underlie Cell Adhesion in Mammalian Cells](https://doi.org/10.3389/fcell.2021.642212)

# 3. Important notes
## 3.1. Nextflow caching and storage capacity
Nextflow creates a lot of temporary files! Even though this is sometimes useful (e.g. when you change a small part of the pipeline, and run the pipeline again with -resume option, it will rerun only the parts that have been affected), it is advised to keep an eye on the directory (`work/`) and clean it from time to time.

## 3.2. Results directory structure
For each dataset processed with this pipelines (and subsequent notebooks) the results directory structure is as follows:

* $dataset
  * clean_panther4march -- family-wise results
    * intermediate
      * dichotomised_genes.csv
      * family_thresholds.csv
      * familyset.csv -- filtered familyset for this specific dataset
    * results
      * family_IC.csv
      * shuffled_family_IC.csv
  * intermediate
    * $dichotomisation_type
      * dichotomised_genes.csv
      * filtered_dichotomised_genes.csv
      * optimal_thresholds.csv
    * chr*_filtered.csv
  * results 
    * $dichotomisation_type
      * shuffled_IC.csv
      * stage1_chr*_IC.csv
  * mathematica -- optional, results were put there manually
    * mathematica_results.csv -- the data obtained from Mathematica
    * fitted_distributions_dichotomised.csv
    * fitted_distributions_thresholds.csv
  * chr_ks.csv -- Kolmogorov-Smirnov test, chromosome-wise
  * chr_stat_test_pvalues*.csv -- permutation tests, chromosome-wise
  * family_stat_test_pvalues.csv -- statistical tests, family-wise


# 4. Software description
## 4.1. Data Preparation
### 4.1.1. Prepare family set (PANTHER): **preprocess_panther.py**
Input: raw PANTHER families file (i.e. PTHR15.0_mouse.tsv)

Output: 
* Intermediate -- preprocessed_panther.csv (various extracted gene ids to be converted into gene symbols)
* Final -- clean_panther.csv (final table)

Description:
Extracts MGI, Ensembl, Uniprot and Entrez IDs. Which are then to be manually converted into gene symbols (e.g. with db2db) and to be again fed into the script. It, then, merges the family IDs with gene symbols. The log of the cleaning that I have performed is below:

* Accessed via ftp://ftp.pantherdb.org/sequence_classifications/current_release/PANTHER_Sequence_Classification_files/PTHR15.0_mouse
* 21393 genes total
* Gene IDs are extracted with regex (Gene Symbols, MGI, Uniprot names, Ensembl IDs and EntrezGene ID)
  * MGI -- 21175 genes
  * Uniprot -- 140 genes
  * Ensembl -- 51 genes
  * Gene Symbols -- 26 genes
  * EntrezGene -- 1 gene
* Family IDs are extracted (by splitting subfamily ids)
* The IDs are converted into gene symbols where possible
  * For MGI -- http://www.informatics.jax.org/batch
  * For the rest -- https://biodbnet-abcc.ncifcrf.gov/db/db2db.php
* Genes with no symbol are removed (172 genes)
* One duplicated gene is renamed (to be aligned with PANTHER) -- MGI:2442135: Arhgef18 -> A430078G23Rik
* One duplicated gene had two distinct families associated:
  Ighd (MGI:96447) -- PTHR23266 and PTHR23411. The first one is chosen, as most of the genes with this prefix are in there
* The rest of the duplicates (11 genes) come from the different accession ids (eg. MGI and Ensembl) and had no discrepancy in family IDs
* The total number of genes left -- 21209

Note: clean_panther.csv is provided in this repository

### 4.1.2. Prepare TPM datasets: **clean_dataset.ipynb**
Input:
* ensembl_id_to_gene_symbol_biomart_export.tsv (from Ensembl biomart)
* *.tsv/*.csv (TPM dataset to be cleaned)

Output: *.csv (cleaned TPM dataset)

Description:
Many datasets that we used had duplicated genes (e.g. Dopaminergic). Therefore the duplicates need to be disabiguated. The process applied:
* Map Ensembl IDs in the raw dataset to the gene symbols from the latest reference (this already fixes the vast majority of the duplicates)
* In case there are still some duplicates (i.e. Ensembl IDs themselves are duplicated), there is no way to differentiate between them (the rows probably come from different transcripts), therefore the rows with the largest expression are kept

## 4.2. Chromosome-wise pipelines
Chromosome-wise pipeline is run by **full_dataset_IC.nf** and **full_dataset_IC_no_BI.nf** (with and without Bimodality Index filtering for TPM and UMI respectively). The final results of the pipeline: dichotomised genes, optimal thresholds for each gene, and IC calculations (both for the original distribution and for the shuffled chromosomes).
Input parameters:
* --gff -- path to the GFF reference file
* --dataset -- path to the *cleaned* TPM/UMI dataset file
* --type -- dichotomisation type (*vrs*, *geomean* or *3max*)
* --shuffle_N -- integer, number of shufflings performed for the permutation test

Examples:

```
time ./full_dataset_IC.nf --dataset `realpath Dopaminergic_TPM_clean.csv` --gff `realpath GRCm38.p6_genomic.gff` --type geomean --shuffle_N 1000 -resume;

time ./full_dataset_IC_no_BI.nf --dataset `realpath Alveolar_cells_Type_I_Merged_Batches.csv` --gff `realpath GRCm38.p6_genomic.gff` --type 3max --shuffle_N 1000 -resume;

```

Note: for the ease of use, all the pipelines are run with **run_workflow.sh** script, you can refer to it for more examples of commands to run

Note: all scripts that are run with pipelines are in **bin/** directory

Note: for all pipelines, the results are saved under the directory corresponding to the name of the dataset file provided.

### 4.2.1. Filter chromosomes: **bin/split_chromosomes.py**
Input: GFF file

Output: \$dataset_name/intermediate/chr*_filtered.csv -- CSV for each filtered chromosome separately

Description:
The script reads the GFF file, splits it into chromosomes, extracts genes, filters them (remove predicted genes, and keep genes with the source BestRefSeq), and outputs results as a series of CSV files.

### 4.2.2. Calculate thresholds and dichotomise: **bin/dichotomise.py**
Input:
* dataset -- cleaned dataset filename
* type -- dichotomisation type (*vrs*, *geomean* or *3max*)

Output:
* \$dataset_name/intermediate/\$type/optimal_thresholds.csv -- optimal threshold found for each gene with the respective dichotomisation method
* \$dataset_name/intermediate/\$type/dichotomised_genes.csv -- table of dichotomised genes

Description:

Calculates *individual* thresholds for each gene. For dichotomisation methods see Methods section in the article, or the corresponding code.

### 4.2.3. Filter genes: **bin/filter_genes.py**
Input:
* dataset -- the set of raw and cleaned genes/expression values
* dichotomised_genes -- the set of dichotomised genes from the previous step
* type -- dichotomisation type
  
Output: \$dataset_name/intermediate/\$type/filtered_dichotomised_genes.csv -- filtered table of dichotomised genes

Description:
Performes BI filtering (for TPM datasets). Other types of filtering can potentially be added.

### 4.2.4. Calculate IC: **bin/calculate_IC.py**
Input: 
* (chr_name, chromosome) -- a tuple of a chromosome name and the filtered chromosome table
* dichotomised_genes
* type -- dichotomisation type

Output:
* (chr_name, \$dataset_name/results/\$type/stage1_chr*_IC.csv) -- the resulting IC calculations (for each chromosome separately)

Description:
The process is run for each chromosome separately (potentially in parallel). The calculations are performed for all overlapping stretches along the chromosome. Stretch sizes 7, 14 and 21.
The output rows record the start and end position in the chromosome (gene number, starting from 0, end -- non-inclusive), stretch size, number of genes in a stretch, observed and PB variances, IC, and mean number of expressed genes per cel.

### 4.2.5. Calculate IC for shuffled chromosomes: **bin/calculate_chromosomes_shuffled_IC.py**
Input:
* dichotomised_genes
* chromosomes -- all chromosome files
* shuffle_N

Output: shuffled_IC.csv (in temporary directory)

Description:
Performes the same calculations as **calculate_IC.py**, but shuffles the order of the genes in chromosomes before that. The process is run shuffle_N times.

### 4.2.6. Collect shuffled IC: **bin/collect_chromosomes_shuffled_IC.py**
Input:
* all shuffled IC files
* type -- dichotomisation type

Output: \$dataset_name/results/\$type/shuffled_IC.csv -- combined shuffled ICs

Description:
When all the shuffled calculations are finished, the temporary files are combined together (adding the column to indicate iterations).


## 4.3. Family-wise pipeline
Family-wise pipeline is run by **family_wise_IC.nf**. The final results of the pipeline: dichotomised genes, optimal thresholds for each family, and IC calculations (both for the original distribution and for the shuffled chromosomes).
Input parameters:
* --familyset -- path to the familyset (e.g. PANTHER)
* --dataset -- path to the *cleaned* TPM/UMI dataset file
* --type -- dichotomisation type (supported: *geomean05*, *fm05*, or *one*)
* --shuffle_N -- integer, number of shufflings performed for the permutation test

Examples:

```
time ./family_wise_IC.nf --dataset `realpath Dopaminergic_TPM_clean.csv` --familyset `realpath clean_panther4march.csv` --type geomean05 --shuffle_N 1000 -resume;

time ./family_wise_IC.nf --dataset `realpath klein.csv` --familyset `realpath clean_panther4march.csv` --type fm05 --shuffle_N 1000 -resume;

```

### 4.3.1. Transform family set: **bin/transform_familyset.py**
Input:
* familyset -- clean PANTHER familyset
* dataset -- clean TPM/UMI dataset

Output:
* familyset.csv -- transformed familyset

Description:
Gets the overlap between the familyset and genes in a dataset, and filters out the families of size less than 5.

### 4.3.2. Dichotomise: **bin/dichotomise_families.py**
Input:
* transformed familyset
* dataset
* type -- dichotomisation type (*geomean05*, *fm05* or *one*)

Output:
* \$dataset_name/\$familyset_name/intermediate/dichotomised_genes.csv
* \$dataset_name/\$familyset_name/intermediate/family_thresholds.csv

Description:
Calculates family-wise thresholds and dichotomises genes accordingly.

### 4.3.3. Calculate IC: **bin/calculate_family_IC.py**
Input:
* transofmed familyset
* dataset
* type

Output: \$dataset_name/\$familyset_name/results/family_IC.csv

Description:
Calculates IC for each family (from dichotomised genes). Reports family ID, family size, number of recorded genes, observed variance, PB variance, IC, and mean number of expressed genes per cell.

### 4.3.4. Permute families: **bin/permute_families.py**
Input:
* transformed familyset
* dataset
* n -- number of shufflings to perform

Output: shuffled_families.csv (temporary)

Description:
The process run N times. Randomly shuffles assignment of genes to gene families.

### 4.3.5. Calculate shuffled IC: **bin/calculate_family_IC.py**
Input: 
* shuffled familyset
* dataset
* type -- dichotomisation type

Output: family_IC.csv (temporary)

Description:
Performs the same calculation and reports the same output as the "Calculate IC" process, but on the shuffled familyset. The process is then run for each shuffled familyset (i.e. N times).

### 4.3.6. Collect Shuffled IC: **bin/collect_shuffled_IC.py**
Input:
* all shuffled_family_IC files

Output: \$dataset_name/\$familyset_name/results/shuffled_family_IC.csv

Description:
Combines all shuffled familyset IC files together, adding the column to indicate iteration.

## 4.4. Chromosome-wise analysis
All of the below is to be run after the pipelines finished calculation, since most of the notebooks rely on the calculation results

### 4.4.1. Statistical tests: **permutation_test.ipynb**
Input:
* dataset name
* dichotomisation type

Description:
Performs Kolmogorov-Smirnov test and permutation test (for each stretch size separately). Plots the results, and saves the results as CSV in \$dataset/ directory

### 4.4.2. Plot statistical tests: **plot_chr_stat_tests.ipynb**
Input:
* results of statistical tests for all datasets

Description:
Combines information on statistical tests for all datasets and plots it (volcano and bar plots). The combined data is then saved in *chr_stats_test_summary.csv*

### 4.4.3. Plot thresholds/dichotomisation: **plot_thresholds.ipynb**
Input:
* dichotomisation information

Description:
Plots histograms marking different thresholds

### 4.4.4. Overlap of tails of chromosomal stretches: **chromosome_overlap_analysis.ipynb**
**Depends on the results of dataset_overlaps.ipynb, run it first**
Input:
* overlap_chromosome_lower_raw.pickle -- from dataset_overlaps.ipynb

Description:
Produces table of overlaps (with verbose description) for chromosomal stretches -- *overlap_chromosome_lower_analysis.xlsx*

## 4.5. Family-wise analysis
### 4.5.1. Statistical tests: **family_wise_analysis.ipynb**
Description:
Performs Kolmogorov-Smirnov test and permutation test (for each family size separately, for family sizes that have at least 30 families). Plots the results, and saves the results as CSV in \$dataset/ directory

### 4.5.2. Overlap of tails of families: **family_overlap_analysis.ipynb**
**Depends on the results of dataset_overlaps.ipynb, run it first**
Input:
* overlap_family_lower_raw.csv -- from dataset_overlaps.ipynb

Description:
Produces table of overlaps (with verbose description) for chromosomal stretches -- *overlap_family_lower_analysis.xlsx*

### 4.5.3. Plot statistical tests: **plot_family_stat_tests.ipynb**
Input:
* results of statistical tests for all datasets

Description:
Combines information on statistical tests for all datasets and plots it (volcano and bar plots). The combined data is then saved in *family_stats_test_summary.csv*

### 4.5.4. Plot chosen families: **plot_family_overlap_heatmap.ipynb**
Plots heatmaps (on/off states) and histograms of number of expressed genes per cell for the chosen families.

## 4.6. Both chromosome- and family-wise analysis
Some notebooks perform analysis on both chromosome- and family-wise calculation results

### 4.6.1. Cell type overlaps: **dataset_overlaps.ipynb**
Calculates lower and upper tail overlaps for both chromosome- and family-wise.

Output:
* overlap_family_upper.xlsx and overlap_family_lower.xlsx
* overlap_family_lower_raw.csv
* overlap_chromosome_upper.xlsx and overlap_chromosome_lower.xlsx
* overlap_chromosome_lower_raw.pickle


## 4.7. Misc
### 4.7.1. Transform Fitted Distributions data (from Mathematica): **fitted_distributions.ipynb**
Transoforms data obtained with Mathematica for further analysis: *\$dataset/mathematica/fitted_distributions_thresholds.csv* and *\$dataset/mathematica/fitted_distributions_dichotomised.csv*

### 4.7.2. Dichotomisation Correlation: **dichotomisation_correlation.ipynb**
Input:
* Thresholds (for all dichotomisation methods and fitted distributions)
* Dichotomised datasets

Description:
Calculates correlation between dichotomisation methods, and plots pair-scatterplots and tables of correlation

### 4.7.3. Thresholds for UMI (no bimodality): **alveolar_exploratory.ipynb**
Plots histograms of marker genes for Alveolar cells, and calculates accuracy of different dichotomisation methods

### 4.7.4. Distribution fitting: **mathematica_code/Dichotomisation.nb**
Input:
* single-cell RNAseq dataset

Description:
Performs distribution fitting for individual genes, and searches for antimodes (if a distribution is bimodal)

### 4.7.5. Permutation test: **mathematics_code/permutation_test.nb**
Input:
* Table of dichotomised genes

Description:
Reimplementation of the shuffling and the permutation test in Mathematica (family-wise only). Does both steps (shuffling and test) and produces the p-values as an output.

## 4.8. Utility functions
### 4.8.1. bin/pipeline_utils/gff_utils.py
* read_gff3(filename) -- reads GFF file into pandas Dataframe
* get_genes_from_gff3(gff3) -- extracts all genes from the dataframe

### 4.8.2. bin/pipeline_utils/ic_utils.py
* get_observed_variance(genes) -- observed variance
* get_pb_variance(genes) -- PB/expected variance
* get_IC(obs_var, pb_var) -- get IC (yes, it just divides one by another, but I kept it for semantics)
* bootstrap_IC(genes, bootstrap_n = 10000) -- can be used for bootstrapping IC value for a set of genes

### 4.8.3. bin/pipeline_utils/nextflow_utils.py
Couple functions that are used by nextflow pipeline specifically

* split_filenames(filenames) -- splits string of form `[filename1, filename2, ...]` into the list of filenames
* map_chromosome_filename(filename) -- get chromosome name/number from the filename (looks for pattern like `chr18`)

### 4.8.4. bin/pipeline_utils/pomegranate_parser.py
Parses Mathematica probability distributions (as text strings) into pomegranate objects (so they can be operated from python)
* to_pomegranate(mathematica) -- converts as described

### 4.8.5. bin/pipeline_utils/vrs_utils.py
* grange(a,b,step=1.2) -- Returns a list between a and b of geometrically progressing series 
* get_vrs(gene, threshold) -- gives VRS for a given gene and threshold
* get_cvrs(gene, threshold) -- gives CVRS for a given gene and threshold

