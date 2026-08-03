import gwaslab as gl
import numexpr as ne
import os
import glob
import logging
from pathlib import Path
import matplotlib.pyplot as plt


########################### GWAS plot (FEV1/FVC) ##############################

mysumstats = gl.Sumstats(
    "/datg/xuxiaopeng/sc_eQTL/07_GWAS/processed_sumstats2/FEV1-FVC.tsv",
            rsid="rsID",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            eaf="EAF",
            p = "P",
            z="Z",
            n="N",
    build="38"
)

width_cm = 18
height_cm = 14

mysumstats.plot_mqq(mode="r",
                    region=(7,99012504,101012504),
                    region_ref = "rs73158411", 
                    region_lead_grid = "False",
                    vcf_path="/share/home/xuxiaopeng/.gwaslab/EAS.ALL.split_norm_af.1kg_30x.hg38.vcf.gz"
                   )

fig = plt.gcf()

fig.set_size_inches(width_cm / 2.54, height_cm / 2.54)

fig.savefig("/datg/xuxiaopeng/sc_eQTL/07_GWAS/coloc_result2/ZKSCAN1_FEV1_FVC.pdf", 
            dpi=300, 
            facecolor="white", 
            bbox_inches="tight")

plt.close()

########################### eQTL plot (AT2a-ZKSCAN1) ##############################

AT2a_ZKSCAN1 = gl.Sumstats(
    "/datg/xuxiaopeng/sc_eQTL/07_GWAS/coloc_result2/locus/AT2a_ENSG00000106261.17_eQTL.txt",
            rsid="variant_id",
            chrom="chrom",
            pos="pos",
            ea="alt",
            nea="ref",
            p = "pval_nominal",
    build="38"
)


width_cm = 18
height_cm = 14

AT2a_ZKSCAN1.plot_mqq(mode="r",
                    region=(7,99012504,101012504),
                    region_ref = "rs73158411", 
                    region_lead_grid = "False",
                    vcf_path="/share/home/xuxiaopeng/.gwaslab/EAS.ALL.split_norm_af.1kg_30x.hg38.vcf.gz"
                   )

fig = plt.gcf()

fig.set_size_inches(width_cm / 2.54, height_cm / 2.54)

fig.savefig("/datg/xuxiaopeng/sc_eQTL/07_GWAS/coloc_result2/ZKSCAN1_AT2a.pdf", 
            dpi=300, 
            facecolor="white", 
            bbox_inches="tight")

plt.close()


########################### GWAS plot (FEV1) ##############################

mysumstats = gl.Sumstats(
    "/datg/xuxiaopeng/sc_eQTL/07_GWAS/processed_sumstats2/FEV1.tsv",
            rsid="rsID",
            chrom="CHR",
            pos="POS",
            ea="EA",
            nea="NEA",
            eaf="EAF",
            p = "P",
            z="Z",
            n="N",
    build="38"
)

width_cm = 18
height_cm = 14

mysumstats.plot_mqq(mode="r",
                    region=(10, 102835672, 104835672),
                    region_ref = "rs10883922", 
                    region_lead_grid = "False",
                    vcf_path="/share/home/xuxiaopeng/.gwaslab/EAS.ALL.split_norm_af.1kg_30x.hg38.vcf.gz"
                   )

fig = plt.gcf()

fig.set_size_inches(width_cm / 2.54, height_cm / 2.54)

fig.savefig("/datg/xuxiaopeng/sc_eQTL/07_GWAS/coloc_result2/STN1_FEV1.pdf", 
            dpi=300, 
            facecolor="white", 
            bbox_inches="tight")

plt.close()

########################### eQTL plot (AT2a-STN1) ##############################

AT2a_STN1 = gl.Sumstats(
    "/datg/xuxiaopeng/sc_eQTL/07_GWAS/coloc_result2/locus/AT2a_ENSG00000107960.12_eQTL.txt",
            rsid="variant_id",
            chrom="chrom",
            pos="pos",
            ea="alt",
            nea="ref",
            p = "pval_nominal",
    build="38"
)

width_cm = 18
height_cm = 14

AT2a_STN1.plot_mqq(mode="r",
                    region=(10, 102835672, 104835672),
                    region_ref = "rs10883922",
                    region_lead_grid = "False",
                    vcf_path="/share/home/xuxiaopeng/.gwaslab/EAS.ALL.split_norm_af.1kg_30x.hg38.vcf.gz"
                   )

fig = plt.gcf()

fig.set_size_inches(width_cm / 2.54, height_cm / 2.54)

fig.savefig("/datg/xuxiaopeng/sc_eQTL/07_GWAS/coloc_result2/STN1_AT2a.pdf", 
            dpi=300, 
            facecolor="white", 
            bbox_inches="tight")

plt.close()
