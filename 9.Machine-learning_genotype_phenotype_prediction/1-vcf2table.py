import os 
import numpy as np
import pandas as pd

from cyvcf2 import VCF

from mpire import WorkerPool

IN_PATH = "/data2/Ricky/risheng/wenzhang/data/1.data/"
OUT_PATH = "/data2/zhoujb/project/hhf/250718/rawData/rawVCFTable/"

def vcf2table(in_file, out_file):
    raw_vcf = VCF(in_file, strict_gt=True, threads=None)

    output_list = []
    output_list.append(["CHROM", "POS", "ID", "REF", "ALT"]+raw_vcf.samples)
    for variant in raw_vcf:
        tmp_out = [variant.CHROM, str(variant.POS), "{}-{}".format(variant.CHROM, variant.POS), variant.REF, ",".join(variant.ALT)]
        for item in variant.genotypes:
            if (item[0] < 0) or (item[1] < 0):
                tmp_out.append("-1")
            else:
                tmp_out.append("{}".format(item[0]+item[1]))
                #if item[0] != item[1]:
                    #tmp_out.append("0")
                #else:
                    #if item[0] == 0:
                        #tmp_out.append("1")
                    #else:
                        #tmp_out.append("-1")

        output_list.append(tmp_out)

    with open(out_file, "w") as out_f:
        for item in output_list:
            print("\t".join(item), file=out_f)

    print("{} DONE".format(out_f))
    return

vcf2table(os.path.join(IN_PATH, "10w/311.new.sample.snp.indel.g95maf05.random100k.vcf.gz"), os.path.join(OUT_PATH, "10w.txt"))
vcf2table(os.path.join(IN_PATH, "150gene/gene150_up2k.311.new.sample.snp.indel.g95maf05.vcf.gz"), os.path.join(OUT_PATH, "150gene.txt"))
vcf2table(os.path.join(IN_PATH, "4gene2k/gene4_up2k.311.new.sample.snp.indel.g95maf05.vcf.gz"), os.path.join(OUT_PATH, "4gene2k.txt"))
vcf2table(os.path.join(IN_PATH, "4gene_updown20k/4gene_updown20k.vcf.gz"), os.path.join(OUT_PATH, "4gene_updown20k.txt"))
vcf2table(os.path.join(IN_PATH, "Fstgene/geneFst_up2k.311.new.sample.snp.indel.g95maf05.vcf.gz"), os.path.join(OUT_PATH, "Fstgene.txt"))

vcf2table(os.path.join(IN_PATH, "GWASvcf/DEC1.MLM.vcf.gz"), os.path.join(OUT_PATH, "GWAS_DEC1.txt"))
vcf2table(os.path.join(IN_PATH, "GWASvcf/DEC2.MLM.vcf.gz"), os.path.join(OUT_PATH, "GWAS_DEC2.txt"))
vcf2table(os.path.join(IN_PATH, "GWASvcf/DEC_BLUE.MLM.vcf.gz"), os.path.join(OUT_PATH, "GWAS_DEC_BLUE.txt"))
vcf2table(os.path.join(IN_PATH, "GWASvcf/PGWC1.MLM.vcf.gz"), os.path.join(OUT_PATH, "GWAS_PGWC1.txt"))
vcf2table(os.path.join(IN_PATH, "GWASvcf/PGWC2.MLM.vcf.gz"), os.path.join(OUT_PATH, "GWAS_PGWC2.txt"))
vcf2table(os.path.join(IN_PATH, "GWASvcf/PGWC_BLUE.MLM.vcf.gz"), os.path.join(OUT_PATH, "GWAS_PGWC_BLUE.txt"))
print("ALL DONE")

