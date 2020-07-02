task abf {
    File zfile
    String pheno
    String prefix = basename(zfile, ".z")
    Float prior_variance
    Float coverage
    String zones
    String docker
    Int cpu
    Int mem

    command <<<
        run_abf.R \
            --z ${zfile} \
            --snp ${prefix}.abf.snp \
            --cred ${prefix}.abf.cred \
            --log ${prefix}.abf.log \
            --prior-variance ${prior_variance} \
            --coverage ${coverage}

    >>>

    output {
        File log = prefix + ".abf.log"
        File snp = prefix + ".abf.snp"
        File cred = prefix + ".abf.cred"

    }

    runtime {

        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 5 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: false
    }
}

task susie {
    Int n_samples
    Float var_y
    Int n_causal_snps
    File zfile
    String pheno
    String prefix = basename(zfile, ".z")
    String zones
    String docker
    Int cpu
    Int mem
    Float min_cs_corr

    command <<<
        #!/usr/bin/env bash

        run_susieR.R \
            --z ${zfile} \
            -n ${n_samples} \
            --L ${n_causal_snps} \
            --var-y ${var_y} \
            --snp ${prefix}.susie.snp \
            --cred ${prefix}.susie.cred \
            --log ${prefix}.susie.log \
            --susie-obj ${prefix}.susie.rds \
            --save-susie-obj \
            --write-alpha \
            --write-single-effect \
            --min-cs-corr ${min_cs_corr}

    >>>

    output {
        File log = prefix + ".susie.log"
        File snp = prefix + ".susie.snp"
        File cred = prefix + ".susie.cred"
        File rds = prefix + ".susie.rds"

    }

    runtime {

        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 5 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: false
    }
}

task combine {
    String pheno
    Int n_causal_snps
    Array[File] abf_snp
    Array[File] abf_cred
    Array[File] susie_snp
    Array[File] susie_cred
    String zones
    String docker
    Int cpu
    Int mem

    command <<<

        cat << "__EOF__" > combine_snp.awk
        BEGIN {
            OFS = "\t"
        }
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                col[$i] = i
            }
            gsub(" ", "\t")
            print "trait", "region", "variant", $0
        }
        FNR == 1 {
            match(FILENAME, /(chr[0-9X]+)\.([0-9]+-[0-9]+)\./, a)
            region = a[1]":"a[2]
        }
        FNR > 1 {
            chrom = substr($col["chromosome"], 4)
            sub(/^0/, "", chrom)
            v = sprintf( \
                "%s:%s:%s:%s", \
                chrom, \
                $col["position"], \
                $col["allele1"], \
                $col["allele2"] \
            )
            gsub(" ", "\t")
            print pheno, region, v, $0 | "sort -V -k2,3"
        }
        __EOF__


        cat << "__EOF__" > combine_cred.awk
        BEGIN {
            OFS = "\t"
        }
        NR == 1 {
            print "trait", "region", $0
        }
        FNR == 1 {
            match(FILENAME, /(chr[0-9X]+)\.([0-9]+-[0-9]+)\./, a)
            region = a[1]":"a[2]
        }
        FNR > 1 {
            print pheno, region, $0
        }
        __EOF__

        # Combine abf .snp files
        awk -f combine_snp.awk -v pheno=${pheno} ${sep=" " abf_snp} | bgzip -c -@ ${cpu} > ${pheno}.ABF.snp.bgz
        tabix -s 5 -b 6 -e 6 -S 1 ${pheno}.ABF.snp.bgz
        # Combine abf .cred files
        awk -f combine_cred.awk -v pheno=${pheno}  ${sep=" " abf_cred} | bgzip -c -@ ${cpu} > ${pheno}.ABF.cred.bgz


        # Combine susie .snp files
        awk -f combine_snp.awk -v pheno=${pheno} ${sep=" " susie_snp} | bgzip -c -@ ${cpu} > ${pheno}.SUSIE.snp.bgz
        tabix -s 5 -b 6 -e 6 -S 1 ${pheno}.SUSIE.snp.bgz
        # Combine susie .cred files
        awk -f combine_cred.awk -v pheno=${pheno}  ${sep=" " susie_cred} | bgzip -c -@ ${cpu} > ${pheno}.SUSIE.cred.bgz

    >>>

    output {

        File out_abf_snp = pheno + ".ABF.snp.bgz"
        File out_abf_snp_tbi = pheno + ".ABF.snp.bgz.tbi"
        File out_abf_cred = pheno + ".ABF.cred.bgz"
        File out_susie_snp = pheno + ".SUSIE.snp.bgz"
        File out_susie_snp_tbi = pheno + ".SUSIE.snp.bgz.tbi"
        File out_susie_cred = pheno + ".SUSIE.cred.bgz"

    }

    runtime {

        docker: "${docker}"
        cpu: "${cpu}"
        memory: "${mem} GB"
        disks: "local-disk 30 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: false
    }
}

workflow finemap {

    String zones
    String docker
    String pheno
    Int n_samples
    Float var_y
    Int n_causal_snps
    Array[File] zfiles

    scatter (zfile in zfiles) {

        call abf {
            input: zones=zones, docker=docker, zfile=zfile, pheno=pheno
        }

        call susie {
            input: zones=zones, docker=docker, zfile=zfile, pheno=pheno,
                n_samples=n_samples, var_y=var_y, n_causal_snps=n_causal_snps,

        }
    }

    call combine {
        input: zones=zones, docker=docker, pheno=pheno, n_causal_snps=n_causal_snps,
            abf_snp=abf.snp, abf_cred=abf.cred, susie_snp=susie.snp, susie_cred=susie.cred
    }
}
