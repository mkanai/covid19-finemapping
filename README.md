# Statistical fine-mapping of the COVID-19 HGI meta-analysis
* Analysis date: June 30, 2020
* Latest meta-analysis release: [Release 3 (June 29, 2020)](https://www.covid19hg.org/results/)

## Summary
We conducted statistical fine-mapping of the meta-analysis results, assuming a single causal variant per locus and a shared causal effect across studies. For each genome-wide significant locus (P < 5e-8; 3 Mb window around the lead variant), we applied approximate Bayes factor (ABF) with a prior variance W = 0.04 ([Wakefield, J. 2009](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.20359)) to estimate posterior inclusion probability (PIP) and 95% credible sets. We observed the variants in 95% CS at the 3p21.31 signal are in very tight LD (r2 > 0.9) across multiple populations, suggesting it might be challenging to disentangle them apart statistically. Although ABF estimates PIP = 0.85 for the lead variant of ANA_C2_V2, we note that this seems suspicious given the LD pattern in a population and is likely due to biases from different phenotyping/genotyping/imputation across studies.

## Results
Results include both numerical tables and locuszoom-like plots.

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
* `lbf`: log Bayes factor
* `prob`: posterior inclusion probability (PIP)

### Locuszoom plots
[Locuszoom](http://locuszoom.org/)-like plots were generated using our custom script. For each fine-mapped locus, we plotted the following panels:
1. Manhattan plot
    * x-axis: genomic position (GRCh38)
    * y-axis: -log10(p)
    * color: r2 value to a lead variant (with maximum PIP) in UKBB Europeans from [the Pan-UKBB project](http://pan.ukbb.broadinstitute.org/).
2. PIP plot
    * x-axis: genomic position (GRCh38)
    * y-axis: PIP
    * color: in 95% CS or not
3. r2 plot
    * x-axis: genomic position (GRCh38)
    * y-axis: r2 value to a lead variant (with maximum PIP) in six populations in UKBB from [the Pan-UKBB project](http://pan.ukbb.broadinstitute.org/).
    * color: Population
    * alpha: variants in 95% CS are plotted with full opacity, whereas those not in CS are dimmed.

## Methods
The whole pipeline will be publicly available at https://github.com/mkanai/covid19-finemapping.

## Known limitations
* Assumption: a single causal variant per locus
    * This is mainly due to the logistical limitations of getting appropriate LD panels, ideally in-sample LD from each cohort, to conduct multi-causal variant fine-mapping.
    * Although we think single-causal variant fine-mapping is sufficient for now given the limited signals in meta-analysis, we plan to apply multi-causal variant fine-mapping when the initiative's analysis plan is ready.
* Assumption: a shared causal effect across studies
    * Previous studies suggest this is a reasonable assumption for most complex traits -- but, we are not entirely sure yet for this study.
    * The column `p_het` helps to identify potentially heterogeneous variants based on Cochran's Q test from meta-analysis.
* We note that, even if the above assumptions are valid, biases due to different phenotyping/genotyping/imputation could still affect single-causal variant fine-mapping, and the results should be interpreted cautiously. For example, in `COVID19_HGI_ANA_C2_V2_20200629`, all the variants in 95% CS show nominal heterogeneity (`p_het` < 0.05; min(`p_het`) = 3.2e-4 for `chr3:45834967:G:GA`), suggesting the observed difference in PIP/marginal p-value is likely due to such biases given the LD pattern in a population.

## Authors
* Masahiro Kanai (mkanai@broadinstitute.org)
* Hilary Finucane (finucane@broadinstitute.org)
