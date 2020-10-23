# Statistical fine-mapping of the COVID-19 HGI meta-analysis
* Analysis date: Oct 22, 2020
* Latest meta-analysis release: [Release 4 (Oct 20, 2020)](https://www.covid19hg.org/results/)

## Summary
We conducted statistical fine-mapping of the meta-analysis results (full genome-wide results from the `leave_23andme` file), assuming a single causal variant per locus and a shared causal effect across studies. For each locus with P < 5e-8 (3 Mb window around the lead variant), we applied approximate Bayes factor (ABF) with a prior variance W = 0.04 ([Wakefield, J. 2009](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.20359)) to estimate posterior inclusion probability (PIP) and 95/99% credible sets. Although ABF estimates several variants with high PIP, **we note that the results should be interpreted cautiously given potential biases from different phenotyping/genotyping/imputation across studies.** In particular, missing variants and heterogeneity among the cohorts could undermine meta-analysis fine-mapping significantly, as illustrated by the discrepancy between a manhattan plot vs. observed LD structure in gnomAD.

## Results
Results include both numerical tables and locuszoom-like plots. Please download them from [the COVID-19 HGI website](https://www.covid19hg.org/).

### Table headers
* `trait`: the whole analysis name including the release date (e.g., `COVID19_HGI_ANA_B2_V2_20200629`)
* `region`: fine-mapped region (`chr:start-end`)
* `variant`: variant id (`chr:position:ref:alt`)
* `rsid`: rsid (to be updated; currently, `chr_position_ref_alt`)
* `chromosome`: chromosome (GRCh38)
* `position`: position (GRCh38)
* `allele1`: reference allele
* `allele2`: alternative allele (effect allele)
* `beta`: marginal beta (log OR) from meta-analysis (`all_inv_var_meta_beta` field)
* `se`: standard error of marginal beta from meta-analysis (`all_inv_var_meta_sebeta` field)
* `p`: p-value from meta-analysis (`all_inv_var_meta_p` field)
* `p_het`: Cochran's Q heterogeneity test p-value from meta-analysis (`all_inv_var_het_p` field)
* `n_studies`: number of studies from meta-analysis (`all_meta_N` field)
* `cs`: 95% credible set (`1` if in CS and `-1` if not)
* `cs_99`: 99% credible set (`1` if in CS and `-1` if not)
* `lbf`: log Bayes factor
* `prob`: posterior inclusion probability (PIP)
* `is_canonical_vep`: Whether the following VEP-annotated consequence is from a canonical transcript.
* `most_severe`: Most severe consequence annotated by VEP
* `gene_most_severe_gene`: Gene symbol that shows most severe consequence
* `gnomad_lead_r2_{pop}`: r2 value in gnomAD to the lead variant in a region

### Locuszoom plots
[Locuszoom](http://locuszoom.org/)-like plots were generated using our custom script. For each fine-mapped locus, we plotted the following panels:
1. Manhattan plot
    * x-axis: genomic position (GRCh38)
    * y-axis: -log10(p)
    * color: (best-guess weighted-average) r2 value to a lead variant (with maximum PIP) which we computed weighted average of r2 value in [gnomAD v2.1.1](https://gnomad.broadinstitute.org/) based on the fraction of ancestral populations in the meta-analysis.
2. PIP plot
    * x-axis: genomic position (GRCh38)
    * y-axis: PIP
    * color: in 99% CS or not
3. r2 plot
    * x-axis: genomic position (GRCh38)
    * y-axis: r2 value to a lead variant (with maximum PIP) in nine populations of [gnomAD v2.1.1](https://gnomad.broadinstitute.org/).
    * color: Population
    * alpha: variants in 99% CS are plotted with full opacity, whereas those not in CS are dimmed.

## Methods
The whole pipeline is publicly available at https://github.com/mkanai/covid19-finemapping.

## Known limitations
* Assumption: a single causal variant per locus
    * This is mainly due to the logistical limitations of getting appropriate LD panels, ideally in-sample LD from each cohort, to conduct multi-causal variant fine-mapping.
    * Although we think single-causal variant fine-mapping is sufficient for now given the limited signals in meta-analysis, we plan to apply multi-causal variant fine-mapping when the initiative's analysis plan is ready.
* Assumption: a shared causal effect across studies
    * Previous studies suggest this is a reasonable assumption for most complex traits -- but, we are not entirely sure yet for this study.
    * The column `p_het` helps to identify potentially heterogeneous variants based on Cochran's Q test from meta-analysis.
* We note that, even if the above assumptions are valid, biases due to different phenotyping/genotyping/imputation could still affect single-causal variant fine-mapping, and the results should be interpreted cautiously. For example, in `COVID19_HGI_ANA_C2_V2_20200629`, all the variants in 99% CS show nominal heterogeneity (`p_het` < 0.05; min(`p_het`) = 3.2e-4 for `chr3:45834967:G:GA`).

## Authors
* Masahiro Kanai (mkanai@broadinstitute.org)
* Hilary Finucane (finucane@broadinstitute.org)
