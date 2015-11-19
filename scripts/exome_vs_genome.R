# Copyright (c) 2015 Genome Research Ltd.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
# of the Software, and to permit persons to whom the Software is furnished to do
# so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

library(mupit)
library(argparse)
library(reshape)
library(ggplot2)
library(Cairo)

get_options <- function() {
    parser = ArgumentParser(description="test the power to detect dominant genes with exome and genome sequencing")
    parser$add_argument("--rates", help="Path to table of mutation rates.",
        default="/nfs/users/nfs_j/jm33/apps/denovonear/results/de_novo_gene_rates.ddd_4k.meta-analysis.txt")
    parser$add_argument("--de-novos", help="Path to DDD de novo dataset.",
        default="/lustre/scratch113/projects/ddd/users/jm33/de_novos.ddd_4k.ddd_only.2015-10-12.txt")
    parser$add_argument("--validations", help="Path to validation results.",
        default="/lustre/scratch113/projects/ddd/users/jm33/de_novos.validation_results.2015-10-12.txt")
    parser$add_argument("--families", help="Path to families PED file.",
        default="/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/family_relationships.txt")
    parser$add_argument("--trios", help="Path to file listing complete trios.",
        default="/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13/trios.txt")
    parser$add_argument("--ddg2p", help="Path to DDG2P file.",
        default="/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2015-04-13/DDG2P/dd_genes_for_clinical_filter")
    parser$add_argument("--iterations", type="integer",
        help="number of iterations to run per condition.", default=100)
    parser$add_argument("--output", help="Path to plot pdf to file.",
        default="exome_vs_genome.pdf")
    
    args = parser$parse_args()
    
    return(args)
}

#' load the table of mutation rates for consequence types per gene
#'
#' @param rates_path path to table of mutation rates per gene.
#'
#' @return dataframe of mutation rates, formatted for running with mupiut
get_rates_dataset <- function(rates_path) {
    rates = read.table(rates_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    # convert from my column names to those used when estimating the gene
    # mutation rates given the cohort size
    rates$hgnc = rates$transcript_id
    rates$mis = rates$missense_rate
    rates$non = rates$nonsense_rate
    rates$splice_site = rates$splice_lof_rate
    rates$syn = rates$synonymous_rate
    rates$frameshift = rates$frameshift_rate
    
    rates = rates[, c("hgnc", "chrom", "length", "mis", "non", "splice_site", "syn", "frameshift")]
    
    return(rates)
}

#' get the probands from trios in the current DDD dataset
#'
#' @param families_path path to family relationships ped file, containing sample
#'        IDs and sex
#' @param trios_path path to table of information on the trios only
#'
#' @return dataframe of proband IDs and sex for the trio-based probands
get_probands <- function(families_path, trios_path) {
    
    families = read.table(families_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    families$is_proband = families$dad_id != "0" | families$mum_id != "0"
    
    # determine the trios with exome data available
    trios = read.table(trios_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    families$is_trio = families$individual_id %in% trios$proband_stable_id
    probands = families[families$is_proband & families$is_trio, ]
    
    probands = probands[, c("individual_id", "sex")]
    
    return(probands)
}

#' get the number of male and female probands in a cohort, or cohort subset
#'
#' @param probands path to table of DDD de novos
#'
#' @return list, with entries for counts of trios with male and female probands
get_trio_counts <- function(probands) {
    # get the number of trios studied in our data for each sex
    sex = table(probands$sex)
    male = sex[["M"]]
    female = sex[["F"]]
    
    return(list("male"=male, "female"=female))
}

#' open de novo mutations from the DDD cohort, removing variants that failed validation
#'
#' @param de_novos_path path to table of DDD de novos
#' @param validations_path path to table of results from validation experiments
#'
#' @return dataframe of de novos in the DDD dataset
get_de_novos <- function(de_novos_path, validations_path) {
    
    diagnosed = NULL
    variants = mupit::get_ddd_de_novos(de_novos_path, diagnosed)
    
    validations = read.table(validations_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    variants = merge(variants, validations, by=c("person_id", "chrom",
        "start_pos", "end_pos", "ref_allele", "alt_allele", "hgnc",
        "consequence"), all.x=TRUE)
    
    # drop out the variants that failed to validate (i.e. were false positives,
    # or inherited)
    variants = variants[!variants$status %in% c("false_positive", "inherited"), ]
    variants$status = NULL
    
    return(variants)
}

#' get the dominant DDG2P genes reaching genomewide significance
#'
#' @param n_trios number of trios to sample from among the probands
#' @param probands dataframe of proband IDs and sex info for DDD probands in trios
#' @param de_novos dataframe of de novos in the DDD dataset
#' @param rates dataframe of mutation rates per gene
#' @param threshold multiple testing corrected threshold for genomewide significance
#' @param genes vector of dominant DDG2P genes
#'
#' @return vector of HGNC symbols for DDG2P genes reaching genomewide significance
get_enrichment_in_sample <- function(n_trios, probands, de_novos, rates, threshold, genes) {
    # Get a sampled subset of probands, then figure out the male and female sex
    # counts in the subset. Restrict the de novos to ones from probands in the
    # subset.
    sampled_ids = try(sample(probands$individual_id, n_trios), silent=TRUE)
    
    # if we try to sample more trios than exists in the current DDD cohort, then
    # return NA, rather than using the full cohort. This will drop out these
    # data points
    if (class(sampled_ids) == "try-error") { return(NULL) }
    
    probands = probands[probands$individual_id %in% sampled_ids, ]
    trios = get_trio_counts(probands)
    de_novos = de_novos[de_novos$person_id %in% sampled_ids, ]
    
    # figure out the enrichment in the
    enrich = mupit::analyse_gene_enrichment(de_novos, trios, rates=rates)
    
    enrich$p_min = apply(enrich[, c("p_lof", "p_func")], 1, min)
    enrich$genomewide = enrich$p_min < threshold
    enrich$dominant = enrich$hgnc %in% genes
    
    dominant = enrich$hgnc[enrich$genomewide & enrich$dominant]
    dominant = dominant[!is.na(dominant)]
    
    return(dominant)
}

#' detect genes reaching genomewide significance, under specific conditions
#'
#' @param probands dataframe of proband IDs and sex info for DDD probands in trios
#' @param de_novos dataframe of de novos in the DDD dataset
#' @param rates dataframe of mutation rates per gene
#' @param threshold multiple testing corrected threshold for genomewide significance
#' @param dominant vector of dominant DDG2P genes
#' @param genomewide vector of genes reaching genomewide significance in
#'        in the complete DDD dataset
#' @param genome_cost reference cost of performing genome sequencing
#' @param budget budget for the current conditions
#' @param relative_cost relative cost of exome sequencing versus genome sequencing
#' @param sensitivity relative sensitivity of genome sequencing for detecting de
#'        novo mutations compared to exome sequencing.
#' @param iterations number of iterations to run for the current conditions
#' @param pb progress bar, to indicate overall progress
#' @param step current step of analysis, for the progress bar
#'
#' @return dataframe of test conditions, numbers of genes reaching genomewide
#'         significance from genome and exome sequencing
run_iterations <- function(probands, de_novos, rates, threshold, dominant,
    genomewide, genome_cost, budget, relative_cost, sensitivity,
    iterations, pb, step) {
    
    # Determine the number of trios that could be sequenced, given the budget.
    # Adjust the number of genome trios upwards by the relative sensitivity of
    # genome sequencing to exome sequencing, to account for the additional de
    # novo mutations that would be detected by genome sequencing.
    n_genome_trios = budget/(genome_cost * 3) * sensitivity
    n_exome_trios = budget/((genome_cost * relative_cost) * 3)
    
    power = data.frame("budget"=numeric(0), "relative_cost"=numeric(0),
        "sensitivity"=numeric(0), "genome"=numeric(0), "exome"=numeric(0))
    n = 0
    while (n < iterations) {
        # We only need to test the exome trios at a sensitivity of 1.0, since
        # variation in sensitity only affects the genome results.
        if (sensitivity == 1.0) {
            exome_genes = get_enrichment_in_sample(n_exome_trios, probands,
                de_novos, rates, threshold, dominant)
        } else {
            exome_genes = NULL
        }
        genome_genes = get_enrichment_in_sample(n_genome_trios, probands,
            de_novos, rates, threshold, dominant)
        
        if (is.null(exome_genes)) { n_exome = NA
        } else { n_exome = sum(exome_genes %in% genomewide) }
        
        if (is.null(genome_genes)) { n_genome = NA
        } else { n_genome = sum(genome_genes %in% genomewide) }
        
        temp = data.frame("budget"=budget, "relative_cost"=relative_cost,
            "sensitivity"=sensitivity, "genome"=n_genome, "exome"=n_exome)
        
        power = rbind(power, temp)
        
        n = n + 1
        step = step + 1
        setTxtProgressBar(pb, step)
    }
    
    return(power)
}

#' simulate power of exome and genome sequencing to detect dominant genes, under various conditions
#'
#' @param probands dataframe of proband IDs and sex info for DDD probands in trios
#' @param de_novos dataframe of de novos in the DDD dataset
#' @param rates dataframe of mutation rates per gene
#' @param threshold multiple testing corrected threshold for genomewide significance
#' @param dominant vector of dominant DDG2P genes
#' @param genomewide vector of genes reaching genomewide significance in
#'        in the complete DDD dataset
#' @param iterations number of iterations to run for the current conditions
#'
#' @return dataframe of power simulations for each condition, containing numbers
#'         of genes reaching genomewide significance from genome and exome sequencing
simulate_power <- function(probands, de_novos, rates, threshold, dominant, genomewide, iterations) {
    budgets = c(1e6, 2e6, 5e6)
    genome_cost = 1000
    exome_relative_cost = seq(0.2, 1.0, 0.2)
    genome_sensitivity = seq(1.0, 1.2, 0.05)
    
    steps = length(budgets) * length(exome_relative_cost) * length(genome_sensitivity) * iterations
    pb = txtProgressBar(min=0, max=steps, initial=NA, style=3)
    
    step = 0
    power = data.frame("budget"=numeric(0), "relative_cost"=numeric(0),
        "sensitivity"=numeric(0), "genome"=numeric(0), "exome"=numeric(0))
    for (budget in budgets) {
        for (relative_cost in exome_relative_cost) {
            for (sensitivity in genome_sensitivity) {
                temp = run_iterations(probands, de_novos, rates, threshold,
                    dominant, genomewide, genome_cost, budget,
                    relative_cost, sensitivity, iterations, pb, step)
                
                power = rbind(power, temp)
                step = step + iterations
            }
        }
    }
    
    return(power)
}

#' plot the results from simulation power of exome and genome sequencing
#'
#' @param power dataframe of power simulations for each condition, containing
#'        numbers of genes reaching genomewide significance from genome and
#'        exome sequencing
#' @param output_path path to save output plot to
plot_power <- function(power, output_path) {
    # reshape the power dataframe, so that we have the mean number of genes
    # reaching genomewide significance for each condition
    power = melt(power, id=c("budget", "relative_cost", "sensitivity"))
    power = cast(power, budget + relative_cost + sensitivity ~ variable, values="value", mean, na.rm=TRUE)
    power = melt(power, id=c("budget", "relative_cost", "sensitivity"))
    power = power[!is.na(power$value), ]
    
    Cairo(output_path, type="pdf", height=15, width=15, units="cm")
    p = ggplot(power, aes(x=relative_cost, y=value, shape=factor(sensitivity), color=variable))
    p = p + facet_grid(. ~ budget)
    p = p + geom_point()
    p = p + geom_line()
    p = p + theme_classic()
    p = p + ggtitle("power of genome and exome sequencing at fixed bugets")
    p = p + ylab("dominant genes at genomewide signficance")
    p = p + xlab("relative cost of exome sequencing to genome sequencing")
    print(p)
    dev.off()
}

main <- function() {
    args = get_options()
    
    rates = get_rates_dataset(args$rates)
    ddg2p = mupit::load_ddg2p(args$ddg2p)
    dominant = sort(unique(ddg2p$gene[ddg2p$mode %in% c("Monoallelic", "X-linked dominant")]))
    
    # set the multiple testing corrected threshold for genomewide significance
    alpha = 0.01
    num_genes = 18500
    num_tests = num_genes * 2
    threshold = alpha/num_tests
    
    probands = get_probands(args$families, args$trios)
    de_novos = get_de_novos(args$de_novos, args$validations)
    genomewide = get_enrichment_in_sample(nrow(probands), probands,
        de_novos, rates, threshold, dominant)
    
    power = simulate_power(probands, de_novos, rates, threshold, dominant,
        genomewide, args$iterations)
    
    plot_power(power, args$output)
}

main()
