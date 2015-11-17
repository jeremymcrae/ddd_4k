# pipeline of processing and analysis scripts for the DDD 4K analyses

# DATE=`date +%Y-%m-%d`
DATE="2015-10-12"
DATAFREEZE="/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13"
DDD1K_DATAFREEZE="/nfs/ddd0/Data/datafreeze/1133trios_20131218/"
USER_DIR="/lustre/scratch113/projects/ddd/users/jm33"
DATA_DIR=${USER_DIR}/"de_novo_data"
FASTA_PATH="/software/ddd/resources/v1.2/hs37d5.fasta"
GENCODE_PATH="/software/ddd/resources/v1.2/gencode.v17.chr_patch_hapl_scaff.annotation.gtf"
DDG2P_PATH="/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2015-04-13/DDG2P/dd_genes_for_clinical_filter"
DE_NOVOS_CALLS_PATH=${DATAFREEZE}/"denovo_gear_trios_extracted_passed_variants_11.05.15.tsv"
FAMILIES_PATH=${DATAFREEZE}/"family_relationships.txt"
TRIOS_PATH=${DATAFREEZE}/"trios.txt"
PHENOTYPES_PATH=${DATAFREEZE}/"phenotypes_and_patient_info.txt"
SAMPLE_IDS_PATH=${DATAFREEZE}/"person_sanger_decipher.txt"
KINSHIP_PATH=${DATAFREEZE}/"kinship_and_pca_trios.txt"
SAMPLE_FAIL_PATH=${DATA_DIR}/"de_novo_sample_fails.txt"
MISSED_INDELS_PATH=${DATA_DIR}/"missed_denovo_indels_datafreeze_2015-04-13.txt"
INDEL_FAILS_PATH=${DATA_DIR}/"de_novo_sample_fails_missed_indels.txt"
DDD_1K_DIAGNOSES=${DDD1K_DATAFREEZE}/"Diagnosis_Summary_1133_20140328.xlsx"
DDD_1K_VALIDATIONS=${DDD1K_DATAFREEZE}/"DNG_Validation_1133trios_20140130.tsv"
DDD_4K_VALIDATIONS=${DATA_DIR}/"de_novos.ddd_4k.validation_results.2015-09-02.xlsx"
LOW_PP_DNM_VALIDATIONS=${DATA_DIR}/"de_novos.ddd_4k.validation_results.low_pp_dnm.2015-10-02.xlsx"

# define paths to put some processed files
LAST_BASE_PATH=${USER_DIR}/"last_base_sites_G.json"
FILTERED_DE_NOVOS_PATH=${USER_DIR}/"de_novos.ddd_4k.ddd_only.${DATE}.txt"
VALIDATIONS_PATH=${USER_DIR}/"de_novos.validation_results.${DATE}.txt"
RATES_PATH=${USER_DIR}/"de_novos.ddd_4k.mutation_rates.${DATE}.txt"
DIAGNOSED_PATH=${USER_DIR}/"ddd_4k.diagnosed.${DATE}.txt"
PHENOTYPES_JSON=${USER_DIR}/"de_novos.ddd_4k.phenotypes_by_proband.json"
MULTISAMPLE_BCF=${USER_DIR}/"ddd_4k.bcftools.bcf"

# define directories for temporary and intermediate files
TEMP_DIR=${USER_DIR}/"temp"
RESULTS_DIR=${USER_DIR}/"results"
AUTOZYGOSITY_DIR=${USER_DIR}/"autozygosity"

# define the paths to the results from the enrichment testing. Some of these
# are input files for later steps
META_WITHOUT_MANHATTAN=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.meta-analysis.manhattan.${DATE}.pdf"
META_WITHOUT_ENRICH=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.meta-analysis.enrichment_results.${DATE}.txt"
META_WITHOUT_CLUSTER=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.meta_analysis.txt"
DDD_WITHOUT_MANHATTAN=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.ddd_only.manhattan.${DATE}.pdf"
DDD_WITHOUT_ENRICH=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.ddd_only.enrichment_results.${DATE}.txt"
DDD_WITHOUT_CLUSTER=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.ddd_only.txt"
DDD_WITHOUT_JSON=${RESULTS_DIR}/"probands_by_gene.without_diagnosed.json"
META_WITH_MANHATTAN=${RESULTS_DIR}/"de_novos.ddd_4k.with_diagnosed.meta-analysis.manhattan.${DATE}.pdf"
META_WITH_ENRICH=${RESULTS_DIR}/"de_novos.ddd_4k.with_diagnosed.meta-analysis.enrichment_results.${DATE}.txt"
META_WITH_CLUSTER=${RESULTS_DIR}/"de_novos.ddd_4k.with_diagnosed.meta_analysis.txt"
DDD_WITH_MANHATTAN=${RESULTS_DIR}/"de_novos.ddd_4k.with_diagnosed.ddd_only.manhattan.${DATE}.pdf"
DDD_WITH_ENRICH=${RESULTS_DIR}/"de_novos.ddd_4k.with_diagnosed.ddd_only.enrichment_results.${DATE}.txt"
DDD_WITH_CLUSTER=${RESULTS_DIR}/"de_novos.ddd_4k.with_diagnosed.ddd_only.txt"
DDD_WITH_JSON=${RESULTS_DIR}/"probands_by_gene.with_diagnosed.json"
ENRICH_WITH_ID=${RESULTS_DIR}/"de_novos.ddd_4k.with_diagnosed.meta_with_ID.enrichment_results.${DATE}.txt"
ENRICH_WITH_ID_AND_AUTISM=${RESULTS_DIR}/"de_novos.ddd_4k.with_diagnosed.meta_with_ID_and_autism.enrichment_results.${DATE}.txt"
ENRICH_WITH_ALL=${RESULTS_DIR}/"de_novos.ddd_4k.with_diagnosed.meta_with_all.enrichment_results.${DATE}.txt"

# define the paths to the results from the de novo clustering
DDD_WITHOUT_CLUSTER_RESULTS=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.ddd_only.clustering_results.txt"
DDD_WITH_CLUSTER_RESULTS=${RESULTS_DIR}/"de_novos.ddd_4k.with_diagnosed.ddd_only.clustering_results.txt"
META_WITHOUT_CLUSTER_RESULTS=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.meta-analysis.clustering_results.txt"
META_WITH_CLUSTER_RESULTS=${RESULTS_DIR}/"de_novos.ddd_4k.with_diagnosed.meta-analysis.clustering_results.txt"

# define the paths to the results from the HPO similarity testing
DDD_WITH_HPOSIMILARITY_RESULTS=${RESULTS_DIR}/"de_novos.ddd_4k.with_diagnosed.ddd_only.phenotype_similarity_results.txt"
DDD_WITHOUT_HPOSIMILARITY_RESULTS=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.ddd_only.phenotype_similarity_results.txt"

# define paths to combined result files
WITH_DIAGNOSED_RESULTS=${RESULTS_DIR}/"de_novos.ddd_4k.with_diagnosed.all.${DATE}.txt"
WITHOUT_DIAGNOSED_RESULTS=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.all.${DATE}.txt"
NOVEL_GENE_VARIANTS=${RESULTS_DIR}/"novel_gene_variants.ddd_4k.${DATE}.txt"
GENES_WITH_DISCREPANT_MECHANISMS=${RESULTS_DIR}/"missing_mechanism_genes.ddd_4k.${DATE}.txt"

# define the paths to candidate CNVs
CANDIDATE_CNVS=${RESULTS_DIR}/"ddd_4k.de_novo_cnvs.${DATE}.txt"
OVERLAPPING_CNVS=${RESULTS_DIR}/"ddd_4k.de_novo_cnvs.overlapping_novel_genes.${DATE}.txt"

################################################################################
# obtain and install required code dependencies
################################################################################

# install necessary python and R packages (ones that are not listed as
# dependencies of the python and R packages on github).
pip install pyfaidx --user
R -e 'install.packages(c("devtools", "argparse"), repos="http://mirrors.ebi.ac.uk/CRAN/")'

# clone all the code repositories for the analyses
git clone https://github.com/jeremymcrae/count_singletons.git
git clone https://github.com/jeremymcrae/denovoFilter.git
git clone https://github.com/jeremymcrae/publishedDeNovos.git
git clone https://github.com/jeremymcrae/mupit.git
git clone https://github.com/jeremymcrae/denovonear.git
git clone https://github.com/jeremymcrae/hpo_similarity.git
git clone https://github.com/jeremymcrae/recessiveStats.git

# install the python packages, a package for filtering de novo calls, a package
# for analysing proximity clustering of de novo mutations, and a package for
# analysing phenotypic similarity between groups of probands
pip install git+git://github.com/jeremymcrae/denovoFilter.git
pip install git+git://github.com/jeremymcrae/denovonear.git
pip install git+git://github.com/jeremymcrae/hpo_similarity.git

# install the R packages, a R package that contains all of the externally
# reported de novo mutations from exome or genome sequencing of trios with
# developmental disorders, and a package to test enrichment of de novo mutations
R -e 'library(devtools);devtools::install_github(jeremymcrae/publishedDeNovos.git)'
R -e 'library(devtools);devtools::install_github(jeremymcrae/mupit.git)'
R -e 'library(devtools);devtools::install_github(jeremymcrae/recessiveStats.git)'

# identify all of the sites where it is possible to have a conserved last base
# of an exon
# runtime: < 10 minutes
python count_singletons/get_last_base_sites.py \
    --fasta ${FASTA_PATH} \
    --gencode ${GENCODE_PATH} \
    --base G \
    --out ${LAST_BASE_PATH}

################################################################################
# filter de novo mutation calls, and prepare validation data
################################################################################

# runtime: < 30 minutes
python denovoFilter/scripts/filter_de_novos.py \
    --de-novos ${DE_NOVOS_CALLS_PATH} \
    --de-novos-indels ${MISSED_INDELS_PATH} \
    --families ${FAMILIES_PATH} \
    --sample-fails ${SAMPLE_FAIL_PATH} \
    --sample-fails-indels ${INDEL_FAILS_PATH} \
    --last-base-sites ${LAST_BASE_PATH} \
    --output ${FILTERED_DE_NOVOS_PATH}

# create a table with all of the validation results in it
# runtime: < 10 minutes
python denovoFilter/scripts/get_validations.py \
    --ddd-1k-validations ${DDD_1K_VALIDATIONS} \
    --ddd-4k-validations ${DDD_4K_VALIDATIONS} \
    --low-pp-dnm ${LOW_PP_DNM_VALIDATIONS} \
    --de-novos ${FILTERED_DE_NOVOS_PATH} \
    --output ${VALIDATIONS_PATH}

################################################################################
# obtain mutation rates for all the genes containing de novo mutations
################################################################################

# prepare a file to create mutation rates from
TEMP_VARS="tmp_variants.txt"
TEMP_GENES="tmp_transcripts.txt"
awk '{ print $7 "\t" $3 "\t" $4 "\t" $9 "\t" $8}' ${FILTERED_DE_NOVOS_PATH} > ${TEMP_VARS}

# identify the transcripts for the genes containing the DDD de novos
# runtime: < 10 hours
python denovonear/scripts/identify_transcripts.py \
    --de-novos ${TEMP_VARS} \
    --minimise-transcripts \
    --out ${TEMP_GENES}

# determine the consequence-specific mutation rates for the genes containing de
# novos
# runtime: < 20 hours, ~ 200 Mb
python denovonear/scripts/construct_mutation_rates.py \
    --genes ${TEMP_GENES} \
    --rates "denovonear/data/forSanger_1KG_mutation_rate_table.txt" \
    --out ${RATES_PATH}

################################################################################
# identify the probands either with diagnoses, or with likely diagnoses
################################################################################

# runtime: < 5 minutes
Rscript mupit/scripts/get_diagnostic_probands.R \
    --ddd-1k-diagnoses ${DDD_1K_DIAGNOSES} \
    --de-novos ${FILTERED_DE_NOVOS_PATH} \
    --low-pp-dnm ${LOW_PP_DNM_VALIDATIONS} \
    --families ${FAMILIES_PATH} \
    --ddg2p ${DDG2P_PATH} \
    --out ${DIAGNOSED_PATH}

################################################################################
# test for enrichment of de novo mutations within genes
################################################################################

# runtime: < 5 minutes
Rscript mupit/scripts/ddd_analysis.R \
    --rates ${RATES_PATH} \
    --de-novos ${FILTERED_DE_NOVOS_PATH} \
    --validations ${VALIDATIONS_PATH} \
    --families ${FAMILIES_PATH} \
    --trios ${TRIOS_PATH} \
    --ddg2p ${DDG2P_PATH} \
    --diagnosed ${DIAGNOSED_PATH} \
    --meta-analysis \
    --meta-subset "intellectual_disability,epilepsy,autism" \
    --out-manhattan ${META_WITHOUT_MANHATTAN} \
    --out-enrichment ${META_WITHOUT_ENRICH} \
    --out-clustering ${META_WITHOUT_CLUSTER}

# runtime: < 5 minutes
Rscript mupit/scripts/ddd_analysis.R \
    --rates ${RATES_PATH} \
    --de-novos ${FILTERED_DE_NOVOS_PATH} \
    --validations ${VALIDATIONS_PATH} \
    --families ${FAMILIES_PATH} \
    --trios ${TRIOS_PATH} \
    --ddg2p ${DDG2P_PATH} \
    --diagnosed ${DIAGNOSED_PATH} \
    --out-manhattan ${DDD_WITHOUT_MANHATTAN} \
    --out-enrichment ${DDD_WITHOUT_ENRICH} \
    --out-clustering ${DDD_WITHOUT_CLUSTER} \
    --out-probands-by-gene ${DDD_WITHOUT_JSON}

# runtime: < 5 minutes
Rscript mupit/scripts/ddd_analysis.R \
    --rates ${RATES_PATH} \
    --de-novos ${FILTERED_DE_NOVOS_PATH} \
    --validations ${VALIDATIONS_PATH} \
    --families ${FAMILIES_PATH} \
    --trios ${TRIOS_PATH} \
    --ddg2p ${DDG2P_PATH} \
    --meta-analysis \
    --meta-subset "intellectual_disability,epilepsy,autism" \
    --out-manhattan ${META_WITH_MANHATTAN} \
    --out-enrichment ${META_WITH_ENRICH} \
    --out-clustering ${META_WITH_CLUSTER}

# runtime: < 5 minutes
Rscript mupit/scripts/ddd_analysis.R \
    --rates ${RATES_PATH} \
    --de-novos ${FILTERED_DE_NOVOS_PATH} \
    --validations ${VALIDATIONS_PATH} \
    --families ${FAMILIES_PATH} \
    --trios ${TRIOS_PATH} \
    --ddg2p ${DDG2P_PATH} \
    --out-manhattan ${DDD_WITH_MANHATTAN} \
    --out-enrichment ${DDD_WITH_ENRICH} \
    --out-clustering ${DDD_WITH_CLUSTER} \
    --out-probands-by-gene ${DDD_WITH_JSON}

################################################################################
# test for enrichment within subsets of the external studies
################################################################################

# runtime: < 5 minutes
Rscript mupit/scripts/ddd_analysis.R \
    --rates ${RATES_PATH} \
    --de-novos ${FILTERED_DE_NOVOS_PATH} \
    --validations ${VALIDATIONS_PATH} \
    --families ${FAMILIES_PATH} \
    --trios ${TRIOS_PATH} \
    --meta-analysis \
    --meta-subset "intellectual_disability" \
    --out-enrichment ${ENRICH_WITH_ID}

# runtime: < 5 minutes
Rscript mupit/scripts/ddd_analysis.R \
    --rates ${RATES_PATH} \
    --de-novos ${FILTERED_DE_NOVOS_PATH} \
    --validations ${VALIDATIONS_PATH} \
    --families ${FAMILIES_PATH} \
    --trios ${TRIOS_PATH} \
    --meta-analysis \
    --meta-subset "intellectual_disability,autism" \
    --out-enrichment ${ENRICH_WITH_ID_AND_AUTISM}

# runtime: < 5 minutes
Rscript mupit/scripts/ddd_analysis.R \
    --rates ${RATES_PATH} \
    --de-novos ${FILTERED_DE_NOVOS_PATH} \
    --validations ${VALIDATIONS_PATH} \
    --families ${FAMILIES_PATH} \
    --trios ${TRIOS_PATH} \
    --meta-analysis \
    --out-enrichment ${ENRICH_WITH_ALL}

python ddd_4k/scripts/check_subset_differences.py \
    --baseline ${DDD_WITH_ENRICH} \
    --modified ${ENRICH_WITH_ID}

python ddd_4k/scripts/check_subset_differences.py \
    --baseline ${DDD_WITH_ENRICH} \
    --modified ${ENRICH_WITH_ID_AND_AUTISM}

python ddd_4k/scripts/check_subset_differences.py \
    --baseline ${DDD_WITH_ENRICH} \
    --modified ${META_WITH_ALL}

python ddd_4k/scripts/check_subset_differences.py \
    --baseline ${DDD_WITH_ENRICH} \
    --modified ${DDD_WITHOUT_ENRICH}

################################################################################
# analyse proximity clustering of de novo mutations
################################################################################

# runtime: < 1 hour (due to being split into job array)
python denovonear/scripts/run_batch.py \
    --script "denovonear/scripts/clustering.py" \
    --temp-dir ${TEMP_DIR} \
    --in ${DDD_WITHOUT_CLUSTER} \
    --rates "denovonear/data/forSanger_1KG_mutation_rate_table.txt" \
    --out ${DDD_WITHOUT_CLUSTER_RESULTS}

# runtime: < 1 hour (due to being split into job array)
python denovonear/scripts/run_batch.py \
    --script "denovonear/scripts/clustering.py" \
    --temp-dir ${TEMP_DIR} \
    --in ${DDD_WITH_CLUSTER} \
    --rates "denovonear/data/forSanger_1KG_mutation_rate_table.txt" \
    --out ${DDD_WITH_CLUSTER_RESULTS}

# runtime: < 1 hour (due to being split into job array)
python denovonear/scripts/run_batch.py \
    --script "denovonear/scripts/clustering.py" \
    --temp-dir ${TEMP_DIR} \
    --in ${META_WITHOUT_CLUSTER} \
    --rates "denovonear/data/forSanger_1KG_mutation_rate_table.txt" \
    --out ${META_WITHOUT_CLUSTER_RESULTS}

# runtime: < 1 hour (due to being split into job array)
python denovonear/scripts/run_batch.py \
    --script "denovonear/scripts/clustering.py" \
    --temp-dir ${TEMP_DIR} \
    --in ${META_WITH_CLUSTER} \
    --rates "denovonear/data/forSanger_1KG_mutation_rate_table.txt" \
    --out ${META_WITH_CLUSTER_RESULTS}

################################################################################
# analyse phenotypic similarity of probands who share de novos in genes
################################################################################

# runtime: < 5 minutes
python hpo_similarity/scripts/prepare_ddd_files.py \
    --phenotypes ${PHENOTYPES_PATH} \
    --sample-ids ${SAMPLE_IDS_PATH} \
    --out ${PHENOTYPES_JSON}

# runtime: < 1 hour (due to being split into job array)
python hpo_similarity/scripts/run_batch.py \
    --script "hpo_similarity/scripts/proband_similarity.py" \
    --temp-dir ${TEMP_DIR} \
    --phenotypes ${PHENOTYPES_JSON} \
    --genes ${DDD_WITHOUT_JSON} \
    --out ${DDD_WITHOUT_HPOSIMILARITY_RESULTS}

# runtime: < 1 hour (due to being split into job array)
python hpo_similarity/scripts/run_batch.py \
    --script "hpo_similarity/scripts/proband_similarity.py" \
    --temp-dir ${TEMP_DIR} \
    --phenotypes ${PHENOTYPES_JSON} \
    --genes ${DDD_WITH_JSON} \
    --out ${DDD_WITH_HPOSIMILARITY_RESULTS}

################################################################################
# merge result files
################################################################################

# runtime: < 5 minutes
Rscript mupit/scripts/combine_all_tests.R \
    --ddg2p ${DDG2P_PATH} \
    --ddd-enrichment ${DDD_WITH_ENRICH} \
    --meta-enrichment ${META_WITH_ENRICH} \
    --ddd-clustering ${DDD_WITH_CLUSTER_RESULTS} \
    --meta-clustering ${META_WITH_CLUSTER_RESULTS} \
    --ddd-phenotype ${DDD_WITH_HPOSIMILARITY_RESULTS} \
    --output ${WITH_DIAGNOSED_RESULTS}

# runtime: < 5 minutes
Rscript mupit/scripts/combine_all_tests.R \
    --ddg2p ${DDG2P_PATH} \
    --ddd-enrichment ${DDD_WITHOUT_ENRICH} \
    --meta-enrichment ${META_WITHOUT_ENRICH} \
    --ddd-clustering ${DDD_WITHOUT_CLUSTER_RESULTS} \
    --meta-clustering ${META_WITHOUT_CLUSTER_RESULTS} \
    --ddd-phenotype ${DDD_WITHOUT_HPOSIMILARITY_RESULTS} \
    --output ${WITHOUT_DIAGNOSED_RESULTS}

################################################################################
# identify variants in candidate novel genes, and prepare reports
################################################################################

python ddd_4k/scripts/get_variants_in_novel_genes.py \
    --de-novos ${FILTERED_DE_NOVOS_PATH} \
    --validations ${VALIDATIONS_PATH} \
    --sanger-ids ${SAMPLE_IDS_PATH} \
    --results ${WITHOUT_DIAGNOSED_RESULTS} \
    --output ${NOVEL_GENE_VARIANTS}

python ddd_4k/scripts/prepare_gene_reports.py \
    --de-novos ${NOVEL_GENE_VARIANTS} \
    --phenotypes ${PHENOTYPES_PATH} \
    --sanger-ids ${SAMPLE_IDS_PATH} \
    --output-dir "gene_reports"

python ddd_4k/scripts/plot_de_novos_in_gene.py \
    --de-novos ${FILTERED_DE_NOVOS_PATH} \
    --external-sites "publishedDeNovos/data-raw/variants.txt.gz" \
    --validations ${VALIDATIONS_PATH} \
    --results ${WITHOUT_DIAGNOSED_RESULTS} \
    --output-dir "gene_reports"

################################################################################
# analyse autozygosity against probability of having a diagnosis
################################################################################

# generate a multi-sample BCF that contains all of the genotype information for
# all of the probands in the DDD
# runtime: < 10 hours (due to being split into multiple jobs)
python recessiveStats/scripts/autozygosity/merge_vcfs.py \
    --vcf-annotate "/software/hgi/pkglocal/vcftools-0.1.11/bin/vcf-annotate" \
    --bcftools "/software/hgi/pkglocal/bcftools-1.2/bin/bcftools" \
    --families ${FAMILIES_PATH} \
    --temp ${TEMP_DIR} \
    --output ${MULTISAMPLE_BCF}

# identify the autozygous regions in all probands
# runtime: < 10 hours (due to being split into multiple jobs)
Rscript mupit/scripts/autozygosity/check_probands_autozygosity.R \
    --script "mupit/scripts/autozygosity/proband_autozygosity.R" \
    --rbinary "/software/R-3.2.2/bin/Rscript" \
    --bcf ${MULTISAMPLE_BCF} \
    --output-folder ${AUTOZYGOSITY_DIR}

# check the likelihood of being diagnostic vs length of autozygous regions
python ddd_4k/scripts/autozygosity_vs_diagnostic.py \
    --autozygosity-dir ${AUTOZYGOSITY_DIR} \
    --consanguinous ${KINSHIP_PATH} \
    --trios ${TRIOS_PATH} \
    --diagnosed ${DIAGNOSED_PATH} \
    --output-groups "ddd_4k/results/autozygosity_vs_diagnosed.groups.pdf"
    --output-regression "ddd_4k/results/autozygosity_vs_diagnosed.regression.pdf"

# run de novo vs phenotypic severity checks
# runtime: < 5 minutes
python ddd_4k/scripts/mutations_by_phenotype.py \
    --de-novos ${FILTERED_DE_NOVOS_PATH} \
    --ddg2p ${DDG2P_PATH} \
    --phenotypes ${PHENOTYPES_PATH} \
    --sanger-ids ${SANGER_IDS} \
    --families ${FAMILIES_PATH} \
    --trios ${TRIOS_PATH} \
    --validations ${VALIDATIONS_PATH} \
    --output-folder "ddd_4k/results"

# identify candidate de novo CNVs in proband VCFs and identify candidate CNVs
# overlapping candidate novel genes
python ddd_4k/get_de_novo_cnvs.py --families ${FAMILIES_PATH} --trios ${TRIOS_PATH} --output ${CANDIDATE_CNVS}
python ddd_4k/get_overlapping_cnvs.py --cnvs ${CANDIDATE_CNVS} --associations ${WITHOUT_DIAGNOSED_RESULTS} --output ${OVERLAPPING_CNVS}

# look into the power of exome sequencing vs genome sequencing
python ddd_4k/scripts/exome_vs_genome.py \
    --output-folder "ddd_4k/results"

# find genes with discrepant mechanisms (i.e. the known genes file lists the
# gene as having a loss-of function mechanism, but we observe missense
# enrichment and clustering, or the known genes file lists a gene as having
# missense variants, but we observe loss-of-fucntion enrichment).
python ddd/scripts/get_genes_with_discrepant_mechanisms.py \
    --known-genes ${KNOWN_GENES} \
    --results ${WITHOUT_DIAGNOSED_RESULTS} \
    --output ${GENES_WITH_DISCREPANT_MECHANISMS}
