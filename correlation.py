import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import statsmodels.stats.multitest as multi


def extract_names(filename: str):
    nameset = []
    with open(filename) as f:
        for index, line in enumerate(f):
            line = line.split("\t")
            if index == 0:
                continue
            nameset.append(np.array(line[0:2]))
    return nameset


# extracting the data from the files
def extract_data(filename: str):
    dataset = []
    with open(filename) as f:
        for index, line in enumerate(f):
            line = line.split("\t")
            if index == 0:
                continue
            dataset.append(np.array(line[2:]).astype(np.float))
    return dataset


data_healthy = extract_data("lusc-rsem-fpkm-tcga_paired.txt")
data_cancerous = extract_data("lusc-rsem-fpkm-tcga-t_paired.txt")
data_names = extract_names("lusc-rsem-fpkm-tcga_paired.txt")


k = 0
for i in data_healthy:
    sum = 0
    for j in i:
        if j == 0:
            sum += 1
    if sum > 25:
        data_healthy.pop(k)
        data_names.pop(k)
        data_cancerous.pop(k)
    k += 1
k = 0
for i in data_cancerous:
    sum = 0
    for j in i:
        if j == 0:
            sum += 1
    if sum > 25:
        data_healthy.pop(k)
        data_names.pop(k)
        data_cancerous.pop(k)
    k += 1


pearson_coeff = []
spearman_coeff = []
# Getting the pearson coeff for each gene
[
    pearson_coeff.append(stats.pearsonr(data_cancerous[i], data_healthy[i])[0])
    for i in range(0, len(data_healthy))
]
pearson_coeff_numpy = np.array(pearson_coeff)
pearson_coeff_numpy[np.isnan(pearson_coeff_numpy)] = 0


# Sort the realtional coefficient array
sorted_coeff = sorted(pearson_coeff_numpy)

# Getting the genes with max and min correlation coeff
max_index = np.argmax(pearson_coeff_numpy)
min_index = np.argmin(pearson_coeff_numpy)
# print(data_healthy[min_index], data_cancerous[min_index])

rel_p_values = []
ind_p_values = []
[
    rel_p_values.append(stats.ttest_rel(data_cancerous[i], data_healthy[i]).pvalue)
    for i in range(0, len(data_healthy))
]
[
    ind_p_values.append(stats.ttest_ind(data_cancerous[i], data_healthy[i]).pvalue)
    for i in range(0, len(data_healthy))
]
# print(rel_p_values)

rel_fdr = multi.multipletests(rel_p_values, method="fdr_bh")
ind_fdr = multi.multipletests(ind_p_values, method="fdr_bh")

# print(rel_p_values)
rel_p_values_h = []
ind_p_values_h = []
for x in rel_p_values:
    if x < 0.05:
        rel_p_values_h.append(True)
        ind_p_values_h.append(True)
    else:
        rel_p_values_h.append(False)
        ind_p_values_h.append(False)

sum = 0
sum2 = 0
for i in range(0, len(rel_fdr[1])):
    if ind_fdr[0][i] != ind_p_values_h[i]:
        sum += 1
    if rel_fdr[0][i] != rel_p_values_h[i]:
        sum2 += 1
