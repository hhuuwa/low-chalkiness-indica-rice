import os 
import numpy as np
import pandas as pd

RAWVCF_PATH = "/data2/zhoujb/project/hhf/250718/rawData/rawVCFTable/"
OUT_PATH = "/data2/zhoujb/project/hhf/250718/rawData/inputData4geneupdown20k/"

phe_data_file = "/data2/Ricky/risheng/wenzhang/data/1.data/phe/311.Chalkiness.phe"
raw_phe_data = pd.read_table(phe_data_file, index_col="Taxa")
phe_col = raw_phe_data.columns

gene_file = pd.read_table(os.path.join(RAWVCF_PATH, "4gene_updown20k.txt"))

# 80
file_list = os.listdir(RAWVCF_PATH)
for file_name in file_list:
    if (file_name == "4gene2k.txt") or (file_name == '4gene_updown20k.txt'):
        tmp_data = pd.read_table(os.path.join(RAWVCF_PATH, file_name))
    else:
        tmp_data = pd.read_table(os.path.join(RAWVCF_PATH, file_name))

        uniq_id = set(tmp_data["ID"].to_list()).difference(gene_file["ID"].to_list())
        tmp_data = tmp_data[tmp_data["ID"].isin(uniq_id)].copy()

        if file_name != "Fstgene.txt":
            tmp_data = tmp_data.sample(n=699, random_state=120).copy()

    tmp_data = tmp_data.set_index("ID")
    tmp_data = tmp_data.drop(columns=['CHROM', 'POS', 'REF', 'ALT'])
    tmp_data = tmp_data.T
    tmp_data = tmp_data.rename(columns={x:"fea_{}".format(x) for x in tmp_data.columns})

    for item in phe_col:
        tmp_phe = raw_phe_data[[item]]
        data_con = pd.concat([tmp_data, tmp_phe], axis=1)
        data_con.to_csv(os.path.join(OUT_PATH, "{}-{}.txt".format(file_name[:-4], item)), sep="\t")

        print("{} DONE".format("{}-{}.txt".format(file_name[:-4], item)))
        del data_con, tmp_phe

    del tmp_data, file_name

print("ALL DONE")
