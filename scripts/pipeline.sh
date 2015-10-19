# pipeline of processing and analysis scripts for the DDD 4K analyses

# DATE=`date +%Y-%m-%d`
DATE="2015-10-12"
DATAFREEZE="/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13"
USER_DIR="/lustre/scratch113/projects/ddd/users/jm33"
FASTA_PATH="/software/ddd/resources/v1.2/hs37d5.fasta"
GENCODE_PATH="/nfs/users/nfs_j/jm33/reference_data/gencode.v19.annotation.gtf.gz"
DE_NOVOS_CALLS_PATH=${DATAFREEZE}/"denovo_gear_trios_extracted_passed_variants_11.05.15.tsv"
MISSED_INDELS_PATH="/nfs/users/nfs_j/jm33/apps/denovoFilter/data/missed_denovo_indels_datafreeze_2015-04-13.txt"
FAMILIES_PATH=${DATAFREEZE}/"family_relationships.txt"
TRIOS_PATH=${DATAFREEZE}/"trios.txt"
PHENOTYPES_PATH=${DATAFREEZE}/"phenotypes_and_patient_info.txt"
SAMPLE_IDS_PATH=${DATAFREEZE}/"person_sanger_decipher.txt"
SAMPLE_FAIL_PATH="/nfs/users/nfs_j/jm33/apps/denovoFilter/data/sample_fails.txt"
INDEL_FAILS_PATH="/nfs/users/nfs_j/jm33/apps/denovoFilter/data/sample_fails_missed_indels.txt"
LAST_BASE_PATH=${USER_DIR}/"last_base_sites_G.json"
FILTERED_DE_NOVOS_PATH=${USER_DIR}/"de_novos.ddd_4k.ddd_only.${DATE}.txt"
DDD_1K_DIAGNOSES="/nfs/ddd0/Data/datafreeze/1133trios_20131218/Diagnosis_Summary_1133_20140328.xlsx"
DDD_1K_VALIDATIONS="/nfs/ddd0/Data/datafreeze/1133trios_20131218/DNG_Validation_1133trios_20140130.tsv"
DDD_4K_VALIDATIONS="/nfs/users/nfs_j/jm33/de_novos.ddd_4k.validation_results.2015-09-02.xlsx"
LOW_PP_DNM_VALIDATIONS="/nfs/users/nfs_j/jm33/de_novos.ddd_4k.validation_results.low_pp_dnm.2015-10-02.xlsx"
DDG2P_PATH="/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2015-04-13/DDG2P/dd_genes_for_clinical_filter"
VALIDATIONS_PATH=${USER_DIR}/"de_novos.validation_results.${DATE}.txt"
RATES_PATH=${USER_DIR}/"de_novos.ddd_4k.mutation_rates.${DATE}.txt"
DIAGNOSED_PATH=${USER_DIR}/"ddd_4k.diagnosed.${DATE}.txt"
PHENOTYPES_JSON=${USER_DIR}"de_novos.ddd_4k.phenotypes_by_proband.json"
TEMP_DIR=${USER_DIR}/"temp"
RESULTS_DIR=${USER_DIR}/"results"

# define the paths to the results from the enrrichment testing. Some of these
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
DDD_WITH_MANHATTAN=${RESULTS_DIR}/"de_novos.ddd_4k.with_diagnosed.meta-analysis.manhattan.${DATE}.pdf"
DDD_WITH_ENRICH=${RESULTS_DIR}/"de_novos.ddd_4k.with_diagnosed.ddd_only.enrichment_results.${DATE}.txt"
DDD_WITH_CLUSTER=${RESULTS_DIR}/"de_novos.ddd_4k.with_diagnosed.ddd_only.txt"
DDD_WITH_JSON=${RESULTS_DIR}/"probands_by_gene.with_diagnosed.json"

# define the paths to the results from the de novo clustering
DDD_WITHOUT_CLUSTER_RESULTS=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.ddd_only.clustering_results.txt"
DDD_WITH_CLUSTER_RESULTS=${RESULTS_DIR}/"de_novos.ddd_4k.with_diagnosed.ddd_only.clustering_results.txt"
META_WITHOUT_CLUSTER_RESULTS=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.meta-analysis.clustering_results.txt"
META_WITH_CLUSTER_RESULTS=${RESULTS_DIR}/"de_novos.ddd_4k.with_diagnosed.meta-analysis.clustering_results.txt"

# define the paths to the results from the HPO similarity testing
DDD_WITH_HPOSIMILARITY_RESULTS=${RESULTS_DIR}/"de_novos.ddd_4k.with_diagnosed.ddd_only.phenotype_similarity_results.txt"
DDD_WITHOUT_HPOSIMILARITY_RESULTS=${RESULTS_DIR}/"de_novos.ddd_4k.without_diagnosed.ddd_only.phenotype_similarity_results.txt"

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

# identify all of the sites where it is possible to have a conserved last base
# of an exon
# runtime: < 10 minutes
python count_singletons/get_last_base_sites.py \
    --fasta ${FASTA_PATH} \
    --gencode ${GENCODE_PATH} \
    --base G \
    --out ${LAST_BASE_PATH}

# run the de novo filtering
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

# obtain mutation rates for all the genes containing de novo mutations
# prepare a file to create mutation rates from
TEMP_VARS="tmp_variants.txt"
TEMP_GENES="tmp_transcripts.txt"
awk '{ print $7 "\t" $3 "\t" $4 "\t" $9 "\t" $8}' ${FILTERED_DE_NOVOS_PATH} > ${TEMP_VARS}

# identify the transcripts for the genes containing the DDD de novos
# runtime: < 5 hours
python denovonear/scripts/identify_transcripts.py \
    --de-novos ${TEMP_VARS} \
    --minimise-transcripts \
    --out ${TEMP_GENES}

# determine the consequence-specific mutation rates for the genes containing de
# novos
# runtime: < 20 hours
python denovonear/scripts/construct_mutation_rates.py \
    --genes ${TEMP_GENES} \
    --rates "denovonear/data/forSanger_1KG_mutation_rate_table.txt" \
    --out ${RATES_PATH}

# identify the probands either with diagnoses, or with likely diagnoses
# runtime: < 5 minutes
Rscript mupit/scripts/get_diagnostic_probands.R \
    --ddd-1k-diagnoses ${DDD_1K_DIAGNOSES} \
    --de-novos ${FILTERED_DE_NOVOS_PATH} \
    --low-pp-dnm ${LOW_PP_DNM_VALIDATIONS} \
    --families ${FAMILIES_PATH} \
    --ddg2p ${DDG2P_PATH} \
    --out ${DIAGNOSED_PATH}

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
    --out-manhattan ${META_WITH_MANHATTAN} \
    --out-enrichment ${META_WITH_ENRICH} \
    --out-clustering ${META_WITH_CLUSTER}

Rscript mupit/scripts/ddd_analysis.R \
    --rates ${RATES_PATH} \
    --de-novos ${FILTERED_DE_NOVOS_PATH} \
    --validations ${VALIDATIONS_PATH} \
    --families ${FAMILIES_PATH} \
    --trios ${TRIOS_PATH} \
    --ddg2p ${DDG2P_PATH} \
    --diagnosed ${DIAGNOSED_PATH} \
    --out-manhattan ${DDD_WITH_MANHATTAN} \
    --out-enrichment ${DDD_WITH_ENRICH} \
    --out-clustering ${DDD_WITH_CLUSTER} \
    --out-probands-by-gene ${DDD_WITH_JSON}

# analyse proximity clustering of de novo mutations
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


# analyse phenotypic similarity of probands who share de novos in genes
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