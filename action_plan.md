## analyses required
- generate forest-like plots for factors influencing de novo burden in DD genes,
  one for binary variables using odds ratios (such as sex and prenatal
  abnormalities) and one for quantitative variables using beta coefficients
  (such as number of HPO terms, paternal age, length of autozygosity).
- modify analysis of phenotypic data inclusion. One plot to compare p-values
  of genomewide genes with and without incorporating HPO similarity, to
  demonstrate that including HPO similarity p-values does not boost any
  additional genes to genomewide significance. The other to compare p-values in
  seizure subset, showing p-values with and without seizure stratification.
- estimate power to detect haploinsufficient genes. Compare mutability
  versus p-values for neurodevelopmental DDG2P genes with loss-of-function
  mechanism. Bin neurodevelopmental loss-of-function genes by clinical
  recognizability and compare number of observed mutations to expected. This
  relationship allows us to estimate the number of loss-of-function genes that
  we don't observe due to ascertainment biases.
- count the genome-wide significant genes with a loss-of-function mechanism.
  These are genomewide genes with 1+ loss-of-function mutations. Given the
  current power to detect loss-of-function genes, estimate total number of
  loss-of-function genes. For example, if we currently have 50% power, then we
  should have found half of the haploinsufficient genes. Compare to number of
  DDG2P loss-of-function genes.
- estimate excess of missense and loss-of-function mutations, given that we
  should expect synonymous to be neutral. Set the synonymous obs/exp ratio to
  1 by identifying a pp_dnm where the ratio is 1. Use this pp_dnm threshold
  for the missense and loss-of-function mutations. This should give us numbers
  for the excess of missense and loss-of-function mutations, which will be one
  estimate of the proportion of probands in our cohort with a pathogenic
  mutation.
- model proportion of pathogenic mutations as loss-of-function versus
  non-loss-of-function with pLI bins. We can examine the models of excess of
  mutations by pLI scores for known haploinsufficient genes, and known
  non-haploinsufficient genes. Our full dataset should be a mixture of those
  two models. We can determine the mixing proportions by assessing a range of
  mixture proportions, and selecting the proportion where the modeled fit best
  fits the observed data. We will also identify the range of viable mixing
  proportions from the variability of model fits near the optimal point.
- estimate excess of mutations in top half of lof and missense pLI bins. This
  will estimate the number of pathogenic loss-of-function mutations, and given
  the proportions we expect to be loss-of-function/non-loss-of-function, we can
  estimate the number of mutations we expect to be non-loss-of-function. This
  number plus the estimated missense with lof mechanism should equal the excess
  of VEP predicted missense mutations. This gives an independent estimate of
  the proportion of the cohort with a pathogenic de novo mutation.
- estimate enrichment factor for the DDD cohort for pathogenic mutations
  compared to the general population. This can be calculated from the obs/exp
  enrichment of loss-of-function SNVs in haploinsufficient known
  neurodevelopmental genes, particularly in the genes less depleted from
  ascertainment bias, i.e. the genes with low recognizability. From there it's
  matter of estimating the proportion of the DDD with any pathogenic de novo,
  which is the proportion with a pathogenic de nov, plus the proportion with a
  pathogenic de novo CNV (~10%), and the proportion missed from our cohort
  (estimated from the missing genes due to recognizability). Once this
  proportion is corrected for the DDD enrichment factor, we have an estimate of
  the prevalence of births with a developmental disorder caused by a pathogenic
  mutation.
- show pathogenic de novo birth prevalence at median paternal age, with younger
  and older estimates from expected additional mutations per year.
- plot distribution of paternal ages in DDD, along with the proportion of
  fathers in different age bins with a child with a diagnostic de novo, which
  is to demonstrate that we do observe a greater prevalence of pathogenic de
  novos in older fathers. Compare DDD paternal age distribution to UK-wide
  paternal age distribution.

I don't understand one page of plots any more, about the percentage of variants
in known DD genes by pLI bin, for loss-of-function and missense mutations. The
sketch expects genes with high pLI scores to have a higher proportion of lof
variants in known DD genes, but that missense percentages will be flat across
different pLI. I'm not sure how this connects to the other analyses. I'll check
back about this once more of the actions are completed.
