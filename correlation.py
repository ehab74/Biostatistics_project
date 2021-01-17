import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import statsmodels.stats.multitest as multi
import pandas as pd

# extracting the Gene name and Gene ID from the files
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


healthy_data = extract_data("lusc-rsem-fpkm-tcga_paired.txt")
cancerous_data = extract_data("lusc-rsem-fpkm-tcga-t_paired.txt")
data_names = extract_names("lusc-rsem-fpkm-tcga_paired.txt")

# Removing the data that has an expression of 0 in more than 50% of the sample
iterator = 0
for row in healthy_data:
    sum = 0
    for column in row:
        if column == 0:
            sum += 1
    if sum > 25:
        healthy_data.pop(iterator)
        data_names.pop(iterator)
        cancerous_data.pop(iterator)
    iterator += 1

iterator = 0
for row in cancerous_data:
    sum = 0
    for column in row:
        if column == 0:
            sum += 1
    if sum > 25:
        healthy_data.pop(iterator)
        data_names.pop(iterator)
        cancerous_data.pop(iterator)
    iterator += 1

# Getting the pearson coeff for each gene
pearson_correlation = []
rel_p_values = []
ind_p_values = []
rel_p_values_alpha = []
ind_p_values_alpha = []


for i in range(0, len(healthy_data)):
    pearson_cc = stats.pearsonr(cancerous_data[i], healthy_data[i])[0]
    pearson_correlation.append(pearson_cc)

    rel_p_value = stats.ttest_rel(cancerous_data[i], healthy_data[i]).pvalue
    ind_p_value = stats.ttest_ind(cancerous_data[i], healthy_data[i]).pvalue
    rel_p_values.append(rel_p_value)
    ind_p_values.append(ind_p_value)

    if rel_p_value < 0.05:
        rel_p_values_alpha.append(True)
    else:
        rel_p_values_alpha.append(False)

    if ind_p_value < 0.05:
        ind_p_values_alpha.append(True)
    else:
        ind_p_values_alpha.append(False)

    if not i:
        max_value = stats.pearsonr(cancerous_data[i], healthy_data[i])[0]
        min_value = max_value
        max_index = 0
        min_index = 0
        continue

    if max_value < pearson_cc:
        max_value = pearson_cc
        max_index = i

    if min_value > pearson_cc:
        min_value = pearson_cc
        min_index = i


pearson_correlation_numpy = np.array(pearson_correlation)
# Export the CC to an excel sheet
df_names = pd.DataFrame(data_names, columns=["Gene name", "Gene ID"])
df_names["Pearson CC"] = pearson_correlation_numpy
filepath = "Correlations.xlsx"
df_names.to_excel(filepath)


# Sort the realtional coefficient array
sorted_coeff = sorted(pearson_correlation_numpy)

# Getting the genes with max and min correlation coeff
# max_index = np.argmax(pearson_correlation_numpy)
# min_index = np.argmin(pearson_correlation_numpy)
# print(healthy_data[min_index], cancerous_data[min_index])

rel_fdr = multi.multipletests(rel_p_values, method="fdr_bh")
ind_fdr = multi.multipletests(ind_p_values, method="fdr_bh")


common_rel_genes = []
common_ind_genes = []
distinct_rel_genes = []
distinct_ind_genes = []


for i in range(0, len(rel_fdr[0])):
    if ind_fdr[0][i] != ind_p_values_alpha[i]:
        distinct_ind_genes.append((data_names[i], ind_fdr[1][i], ind_p_values[i]))
    else:
        common_ind_genes.append((data_names[i], ind_fdr[1][i], ind_p_values[i]))

    if rel_fdr[0][i] != rel_p_values_alpha[i]:
        distinct_rel_genes.append((data_names[i], rel_fdr[1][i], rel_p_values[i]))
    else:
        common_rel_genes.append((data_names[i], rel_fdr[1][i], rel_p_values[i]))

#print(len(distinct_ind_genes),len(distinct_rel_genes))
print(len(distinct_ind_genes),len(distinct_rel_genes))
#Excel sheet for common relative genes
df_com_rel = pd.DataFrame(common_rel_genes, columns=['Gene Name,ID','After FDR','Before FDR']) 
filepath = 'common_rel_genes.xlsx'
df_com_rel.to_excel(filepath, index=False)
#Excel sheet for common independent genes
df_com_ind = pd.DataFrame(common_ind_genes, columns=['Gene Name,ID','After FDR','Before FDR']) 
filepath = 'common_ind_genes.xlsx'
df_com_ind.to_excel(filepath, index=False)
#Excel sheet for distinct relative genes
df_dis_rel = pd.DataFrame(distinct_rel_genes, columns=['Gene Name,ID','After FDR','Before FDR']) 
filepath = 'distinct_rel_genes.xlsx'
df_dis_rel.to_excel(filepath, index=False)
#Excel sheet for distinct indpendent genes
df_dis_ind = pd.DataFrame(distinct_ind_genes, columns=['Gene Name,ID','After FDR','Before FDR']) 
filepath = 'distinct_ind_genes.xlsx'
df_dis_ind.to_excel(filepath, index=False)
# Plotting

fig, (ax1, ax2) = plt.subplots(2)
fig.suptitle(
    "Gene expression for the highest and lowest correlation coefficient genes."
)
ax1.set_title(f"Gene with Maximum Correlation Coefficient{data_names[max_index]}")
ax2.set_title(f"Gene with Minimum Correlation Coefficient{data_names[min_index]}")
ax1.scatter(healthy_data[max_index], cancerous_data[max_index])
ax1.set(xlabel="healthy_sample", ylabel="cancerous_sample")
ax1.plot(
    np.unique(healthy_data[max_index]),
    np.poly1d(np.polyfit(healthy_data[max_index], cancerous_data[max_index], 1))(
        np.unique(healthy_data[max_index])
    ),
    color="blue",
)
ax2.scatter(healthy_data[min_index], cancerous_data[min_index])
ax2.set(xlabel="healthy_sample", ylabel="cancerous_sample")
ax2.plot(
    np.unique(healthy_data[min_index]),
    np.poly1d(np.polyfit(healthy_data[min_index], cancerous_data[min_index], 1))(
        np.unique(healthy_data[min_index])
    ),
    color="blue",
)
ax1.legend(["relation after fitting", "data"])
ax2.legend(["relation after fitting", "data"])
plt.show()
