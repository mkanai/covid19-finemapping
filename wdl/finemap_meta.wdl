import "finemap_meta_sub.wdl" as sub

task preprocess {
    String pheno
    String sumstats_pattern
    File sumstats = sub(sumstats_pattern,"\\{PHENO\\}",pheno)
    File metadata
    String zones
    String docker
    Int cpu
    Int mem
    Boolean scale_se_by_pval
    Boolean x_chromosome
    String rsid_col
    String chromosome_col
    String position_col
    String allele1_col
    String allele2_col
    String freq_col
    String beta_col
    String se_col
    String p_col
    String delimiter
    # can be helpful if adding finemapping with relaxed threshold after more stringent has already ben run.
    # does not include regions with lead snp < this
    Float p_threshold
    Float? minimum_pval

    command <<<

        zcat ${sumstats} | awk '
        BEGIN {
            FS = "\t"
            OFS = "\t"
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                col[$i] = i
            }
            print "${chromosome_col}", "${position_col}", "${allele1_col}", "${allele2_col}", "freq", "${beta_col}", "${se_col}", "${p_col}", "p_het", "n_studies"
        }
        NR > 1 {
            freq = 0.1 # for now
            print $col["${chromosome_col}"], $col["${position_col}"], $col["${allele1_col}"], $col["${allele2_col}"], freq, $col["${beta_col}"], $col["${se_col}"], $col["${p_col}"], $col["all_inv_var_het_p"], $col["all_meta_N"]
        }
        ' > ${pheno}.sumstats.txt

        awk -v pheno=${pheno} '
        BEGIN {
            FS = "\t"
            OFS = "\t"
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                col[$i] = i
            }
        }
        NR > 1 && $col["analysis"] == pheno {
            phi = $col["n_cases"] / $col["n_samples"]
            var_y = phi*(1-phi)
            print $col["n_samples"] > "n_samples.txt"
            print var_y > "var_y.txt"
            exit 1
        }
        ' ${metadata}

        make_finemap_inputs.py \
            --sumstats ${pheno}.sumstats.txt \
            --rsid-col "${rsid_col}" \
            --chromosome-col "${chromosome_col}" \
            --position-col "${position_col}" \
            --allele1-col "${allele1_col}" \
            --allele2-col "${allele2_col}" \
            --freq-col "${freq_col}" \
            --beta-col "${beta_col}" \
            --se-col "${se_col}" \
            --p-col "${p_col}" \
            --extra-cols "p_het" "n_studies" \
            --delimiter "${delimiter}" \
            --set-variant-id \
            --grch38 \
            --exclude-MHC \
            --no-upload \
            --prefix ${pheno} \
            --out ${pheno} \
            --wdl \
            ${true='--scale-se-by-pval ' false=' ' scale_se_by_pval} \
            ${true='--x-chromosome' false=' ' x_chromosome} \
            --p-threshold ${p_threshold} \
            ${true='--min-p-threshold ' false='' defined(minimum_pval)}${minimum_pval}

        res=`cat ${pheno}_had_results`

        if [ "$res" == "False" ]
        then
            touch ${pheno}".z"
            touch ${pheno}".lead_snps.txt"
            touch ${pheno}".lead_snps.txt"
            touch ${pheno}".bed"
        else
            # remove dummy maf field for now
            for z in ${pheno}.*.z
            do
                cp $z $z.old
                cut -d' ' -f6 --complement $z.old > $z
                rm $z.old
            done
        fi

    >>>

    output {

        Array[File] zfiles = glob("*.z")
        File leadsnps = pheno + ".lead_snps.txt"
        File bed = pheno + ".bed"
        File log = pheno + ".log"
        Boolean had_results = read_boolean("${pheno}_had_results")
        Int n_samples = read_int("n_samples.txt")
        Float var_y = read_float("var_y.txt")
    }

    runtime {
        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 20 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: false
    }
}

workflow finemap_meta {

    String zones
    String docker
    String sumstats_pattern
    File phenolistfile

    Array[String] phenos = read_lines(phenolistfile)

    scatter (pheno in phenos) {

        call preprocess {
            input: zones=zones, docker=docker, pheno=pheno, sumstats_pattern=sumstats_pattern
        }

        if( preprocess.had_results) {
            call sub.finemap {
                input: zones=zones, docker=docker, pheno=pheno, zfiles=preprocess.zfiles, pheno=pheno,
                    n_samples=preprocess.n_samples, var_y=preprocess.var_y
            }
        }

    }
}
