# coding: utf-8
import argparse
import numpy as np
import hail as hl
from gnomad.utils.vep import process_consequences
from hail.linalg import BlockMatrix
from hail.utils import new_temp_file


POPS = ['afr', 'amr', 'asj', 'eas', 'est', 'fin', 'nfe', 'nwe', 'seu']
rg37 = hl.get_reference('GRCh37')
rg38 = hl.get_reference('GRCh38')
rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)

# cf: https://github.com/astheeggeggs/BipEx/blob/master/scripts_Dalio/QC_Dalio/07_annotate_variants.py
# https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/0.1/08.categorize_vep_consequences.py
ptv = hl.set(
    ["transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant"])
missense = hl.set([
    "stop_lost", "start_lost", "transcript_amplification", "inframe_insertion", "inframe_deletion", "missense_variant",
    "protein_altering_variant", "splice_region_variant"
])
synonymous = hl.set(["incomplete_terminal_codon_variant", "stop_retained_variant", "synonymous_variant"])
non_coding = hl.set([
    "coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant", "intron_variant", "NMD_transcript_variant", "non_coding_transcript_variant",
    "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification",
    "TF_binding_site_variant", "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation",
    "regulatory_region_variant", "feature_truncation", "intergenic_variant"
])


def annotate_consequence_category(csq_expr, annot_location='consequence_category'):
    annot_expr = {annot_location: hl.case()
                                    .when(ptv.contains(csq_expr), "ptv")
                                    .when(missense.contains(csq_expr), "missense")
                                    .when(synonymous.contains(csq_expr), "synonymous")
                                    .when(non_coding.contains(csq_expr), "non_coding")
                                    .or_missing()}
    return annot_expr


def get_diag_mat(diag_vec: BlockMatrix):
    x = diag_vec.T.to_numpy()
    diag_mat = np.identity(len(x)) * np.outer(np.ones(len(x)), x)
    return BlockMatrix.from_numpy(diag_mat)


def main(args):
    ht_snp = hl.import_table(args.snp, impute=True)
    ht_snp = ht_snp.annotate(variant=hl.delimit(
        [ht_snp.chromosome, hl.str(ht_snp.position), ht_snp.allele1, ht_snp.allele2], delimiter=':'))
    ht_snp = ht_snp.annotate(**hl.parse_variant(ht_snp.variant, reference_genome='GRCh38'))
    ht_snp = ht_snp.key_by('locus', 'alleles')
    ht_snp = ht_snp.add_index('idx_snp')
    ht_snp = ht_snp.checkpoint(new_temp_file())

    # annotate vep
    gnomad = hl.read_table('gs://gnomad-public-requester-pays/release/3.0/ht/genomes/gnomad.genomes.r3.0.sites.ht')
    ht_snp = ht_snp.join(gnomad.select('vep'), how='left')
    ht_snp = process_consequences(ht_snp)

    # extract most severe
    ht_snp = ht_snp.annotate(vep=(hl.case().when(hl.is_defined(ht_snp.vep.worst_csq_for_variant_canonical),
                                                 ht_snp.vep.worst_csq_for_variant_canonical).when(
                                                     hl.is_defined(ht_snp.vep.worst_csq_for_variant),
                                                     ht_snp.vep.worst_csq_for_variant).or_missing()),
                             is_canonical_vep=hl.is_defined(ht_snp.vep.worst_csq_for_variant_canonical))
    ht_snp = ht_snp.annotate(most_severe=hl.if_else(hl.is_defined(ht_snp.vep), ht_snp.vep.most_severe_consequence,
                                                    'intergenic_variant'),
                             gene_most_severe=ht_snp.vep.gene_symbol)
    ht_snp = ht_snp.select_globals()
    ht_snp = ht_snp.drop('vep')
    ht_snp = ht_snp.annotate(**annotate_consequence_category(ht_snp.most_severe))
    ht_snp = ht_snp.checkpoint(new_temp_file())

    df = ht_snp.key_by().drop('locus', 'alleles', 'variant', 'idx_snp').to_pandas()

    # annotate LD
    for pop in POPS:
        ht = hl.read_table(
            f'gs://gnomad-public-requester-pays/release/2.1.1/ld/gnomad.genomes.r2.1.1.{pop}.common.adj.ld.variant_indices.ht')
        ht = ht.annotate(locus_hg38=hl.liftover(ht.locus, 'GRCh38'))
        ht = ht.filter(hl.is_defined(ht.locus_hg38))
        ht = ht.key_by('locus_hg38', 'alleles').drop('locus')
        ht = ht_snp.join(ht, 'inner')
        ht = ht.checkpoint(new_temp_file())

        lead_idx = ht.order_by(hl.desc(ht.prob)).head(1).idx.collect()
        idx = ht.idx.collect()
        bm = BlockMatrix.read(f'gs://gnomad-public-requester-pays/release/2.1.1/ld/gnomad.genomes.r2.1.1.{pop}.common.ld.bm')
        bm = bm.filter(idx, idx)
        # re-densify triangluar matrix
        bm = bm + bm.T - get_diag_mat(bm.diagonal())
        bm = bm.filter_rows(np.where(np.array(idx) == lead_idx[0])[0].tolist()) ** 2

        idx_snp = ht.idx_snp.collect()
        r2 = bm.to_numpy()[0]
        df[f'gnomad_lead_r2_{pop}'] = np.nan
        df[f'gnomad_lead_r2_{pop}'].iloc[idx_snp] = r2

    if args.out.startswith('gs://'):
        fopen = hl.hadoop_open
    else:
        fopen = open

    with fopen(args.out, 'w') as f:
        df.to_csv(f, sep='\t', na_rep='NA', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--snp', type=str, required=True, help='Input snp file from fine-mapping')
    parser.add_argument('--out', type=str, required=True, help='Output path')
    parser.add_argument('--delimiter', type=str, default=' ', help='Delimiter for output ld matrix')

    args = parser.parse_args()

    main(args)
