import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

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
pearson_coeff = []
[pearson_coeff.append(stats.pearsonr(data_cancerous[i],data_healthy[i])[0]) for i in range(0,len(data_healthy))]
pearson_coeff_numpy = np.array(pearson_coeff)
pearson_coeff_numpy[np.isnan(pearson_coeff_numpy)] = 0

#Sort the realtional coefficient array
sorted_coeff = sorted(pearson_coeff_numpy)

#Getting the genes with max and min correlation coeff
max_index = np.argmax(pearson_coeff_numpy)
min_index = np.argmin(pearson_coeff_numpy)


#Plotting
fig, (ax1, ax2) = plt.subplots(2)
fig.suptitle('Gene expression for the highest and lowest correlation coefficient genes.')
ax1.set_title('Maximum Correlation Coefficient')
ax2.set_title('Minimum Correlation Coefficient')
ax1.plot(data_healthy[max_index])
ax1.plot(data_cancerous[max_index])
ax2.plot(data_healthy[min_index])
ax2.plot(data_cancerous[min_index])
ax1.legend(['data_healthy','data_cancerous'])
ax2.legend(['data_healthy','data_cancerous'])
plt.show()