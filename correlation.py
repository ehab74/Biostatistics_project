import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import pandas as pd
#extracting the Gene name and Gene ID from the files
def extract_names(filename: str):
    nameset = []
    with open(filename) as f:
        for index,line in enumerate(f):
            line=line.split('\t')
            if index == 0:
                continue
            nameset.append(np.array(line[0:2]))
    return nameset

#extracting the data from the files
def extract_data(filename: str):
    dataset = []
    with open(filename) as f:
        for index,line in enumerate(f):
            line=line.split('\t')
            if index == 0:
                continue
            dataset.append(np.array(line[2:]).astype(np.float))
    return dataset

data_healthy=extract_data('lusc-rsem-fpkm-tcga_paired.txt')
data_cancerous=extract_data('lusc-rsem-fpkm-tcga-t_paired.txt')
data_names=extract_names('lusc-rsem-fpkm-tcga_paired.txt')
#Removing the data that has an expression of 0 in more than 50% of the sample
k = 0
for i in  data_healthy:
    sum=0
    for j in i:
        if j == 0:
            sum += 1
    if sum > 25:
        data_healthy.pop(k)
        data_names.pop(k)
        data_cancerous.pop(k)
    k += 1
k = 0
for i in  data_cancerous:
    sum=0
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

#Getting the pearson coeff for each gene
[pearson_coeff.append(stats.pearsonr(data_cancerous[i],data_healthy[i])[0]) for i in range(0,len(data_healthy))]
pearson_coeff_numpy = np.array(pearson_coeff)
#Export the CC to an excel sheet
df_names = pd.DataFrame (data_names, columns=['Gene name', 'Gene ID'])
df_names['Pearson CC'] = pearson_coeff_numpy
filepath = 'Correlations.xlsx'
df_names.to_excel(filepath)

#getting spearman coeff
# [spearman_coeff.append(stats.spearmanr(data_cancerous[i],data_healthy[i])[0]) for i in range(0,len(data_healthy))]
# spearman_coeff_numpy = np.array(spearman_coeff)
# spearman_coeff_numpy[np.isnan(spearman_coeff_numpy)] = 0
# sorted_spearman_coeff = sorted(spearman_coeff_numpy)
# spearman_max_index = np.argmax(spearman_coeff_numpy)
# spearman_min_index = np.argmin(spearman_coeff_numpy)
# print(data_healthy[spearman_max_index], data_cancerous[spearman_max_index])

#Sort the realtional coefficient array
sorted_coeff = sorted(pearson_coeff_numpy)

#Getting the genes with max and min correlation coeff
max_index = np.argmax(pearson_coeff_numpy)
min_index = np.argmin(pearson_coeff_numpy)
# print(data_healthy[min_index], data_cancerous[min_index])

#Plotting
fig, (ax1, ax2) = plt.subplots(2)
fig.suptitle('Gene expression for the highest and lowest correlation coefficient genes.')
ax1.set_title(f"Gene with Maximum Correlation Coefficient{data_names[max_index]}")
ax2.set_title(f"Gene with Minimum Correlation Coefficient{data_names[min_index]}")
ax1.scatter(data_healthy[max_index],data_cancerous[max_index])
ax1.set(xlabel='data_healthy', ylabel='data_cancerous')
ax1.plot(np.unique(data_healthy[max_index]), np.poly1d(np.polyfit(data_healthy[max_index], data_cancerous[max_index], 1))(np.unique(data_healthy[max_index])), color='blue')
ax2.scatter(data_healthy[min_index],data_cancerous[min_index])
ax2.plot(np.unique(data_healthy[min_index]), np.poly1d(np.polyfit(data_healthy[min_index], data_cancerous[min_index], 1))(np.unique(data_healthy[min_index])), color='blue')

ax1.legend(['relation after fitting', 'data'])
ax2.legend(['relation after fitting', 'data'])
plt.show()