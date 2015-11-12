# pipeline of processing and analysis scripts for the DDD 4K analyses. Most of
# the steps here require the setup steps of the normal pipeline.sh script to be
# run first.

# DATE=`date +%Y-%m-%d`
DATE="2015-10-12"
DATAFREEZE="/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13"
DDD1K_DATAFREEZE="/nfs/ddd0/Data/datafreeze/1133trios_20131218/"
USER_DIR="/lustre/scratch113/projects/ddd/users/jm33"
DATA_DIR=${USER_DIR}/"de_novo_data"
DDG2P_PATH="/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2015-04-13/DDG2P/dd_genes_for_clinical_filter"
FAMILIES_PATH=${DATAFREEZE}/"family_relationships.txt"
TRIOS_PATH=${DATAFREEZE}/"trios.txt"
PHENOTYPES_PATH=${DATAFREEZE}/"phenotypes_and_patient_info.txt"
SAMPLE_IDS_PATH=${DATAFREEZE}/"person_sanger_decipher.txt"
DDD_1K_DIAGNOSES=${DDD1K_DATAFREEZE}/"Diagnosis_Summary_1133_20140328.xlsx"
LOW_PP_DNM_VALIDATIONS=${DATA_DIR}/"de_novos.ddd_4k.validation_results.low_pp_dnm.2015-10-02.xlsx"

# define paths to put some processed files
LAST_BASE_PATH=${USER_DIR}/"last_base_sites_G.json"
ALL_FILTERED_DE_NOVOS_PATH=${USER_DIR}/"de_novos.ddd_4k.ddd_only.${DATE}.txt"
FILTERED_DE_NOVOS_PATH=${USER_DIR}/"de_novos.ddd_4k.ddd_only.${DATE}.high_quality.txt"
VALIDATIONS_PATH=${USER_DIR}/"de_novos.validation_results.${DATE}.txt"
RATES_PATH=${USER_DIR}/"de_novos.ddd_4k.mutation_rates.${DATE}.txt"
DIAGNOSED_PATH=${USER_DIR}/"ddd_4k.diagnosed.${DATE}.high_quality.txt"
PHENOTYPES_JSON=${USER_DIR}/"de_novos.ddd_4k.phenotypes_by_proband.json"

# define directories for temporary and intermediate files
TEMP_DIR=${USER_DIR}/"temp"
RESULTS_DIR=${USER_DIR}/"results_high_quality"

# define the paths to the results from the enrichment testing. Some of these
# are input files for later steps
META_WITHOUT_MANHATTAN=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.meta-analysis.manhattan.${DATE}.pdf"
META_WITHOUT_ENRICH=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.meta-analysis.enrichment_results.${DATE}.txt"
META_WITHOUT_CLUSTER=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.meta_analysis.txt"
DDD_WITHOUT_MANHATTAN=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.ddd_only.manhattan.${DATE}.pdf"
DDD_WITHOUT_ENRICH=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.ddd_only.enrichment_results.${DATE}.txt"
DDD_WITHOUT_CLUSTER=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.ddd_only.txt"
DDD_WITHOUT_JSON=${RESULTS_DIR}/"probands_by_gene.without_diagnosed.json"

# define the paths to the results from the de novo clustering
DDD_WITHOUT_CLUSTER_RESULTS=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.ddd_only.clustering_results.txt"
META_WITHOUT_CLUSTER_RESULTS=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.meta-analysis.clustering_results.txt"

# define the paths to the results from the HPO similarity testing
DDD_WITHOUT_HPOSIMILARITY_RESULTS=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.ddd_only.phenotype_similarity_results.txt"

# define paths to combined result files
WITHOUT_DIAGNOSED_RESULTS=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.all.${DATE}.txt"
NOVEL_GENE_VARIANTS=${RESULTS_DIR}/"novel_gene_variants.ddd_4k.${DATE}.txt"
GENES_WITH_DISCREPANT_MECHANISMS=${RESULTS_DIR}/"missing_mechanism_genes.ddd_4k.${DATE}.txt"

################################################################################
# identify the high quality candidate de novos
################################################################################

python ddd_4k/scripts/get_high_confidence_de_novos.py \
    --de-novos ${ALL_FILTERED_DE_NOVOS_PATH} \
    --validations ${VALIDATIONS_PATH} \
    --output ${FILTERED_DE_NOVOS_PATH}

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
    --in ${META_WITHOUT_CLUSTER} \
    --rates "denovonear/data/forSanger_1KG_mutation_rate_table.txt" \
    --out ${META_WITHOUT_CLUSTER_RESULTS}

################################################################################
# analyse phenotypic similarity of probands who share de novos in genes
################################################################################

# runtime: < 1 hour (due to being split into job array)
python hpo_similarity/scripts/run_batch.py \
    --script "hpo_similarity/scripts/proband_similarity.py" \
    --temp-dir ${TEMP_DIR} \
    --phenotypes ${PHENOTYPES_JSON} \
    --genes ${DDD_WITHOUT_JSON} \
    --out ${DDD_WITHOUT_HPOSIMILARITY_RESULTS}

################################################################################
# merge result files
################################################################################

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
    --variants ${NOVEL_GENE_VARIANTS} \
    --phenotypes ${PHENOTYPES_PATH} \
    --sanger-ids ${SAMPLE_IDS_PATH} \
    --output "gene_reports"
