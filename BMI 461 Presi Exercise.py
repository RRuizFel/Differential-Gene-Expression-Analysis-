import os                           # Get access to operating system
import pandas as pd                 # For data frame
import numpy as np                  # Calculations
import matplotlib.pyplot as plt     # Graph visualization

# Empty data frame to store data
merged_df = pd.DataFrame()

# List of working files, using pathnames
paths = ['/Users/robertoruizfelix/Downloads/Normal Samples (tumor)/TCGA-CZ-4863.tsv',
         '/Users/robertoruizfelix/Downloads/Normal Samples (tumor)/TCGA-CW-5585.tsv',
         '/Users/robertoruizfelix/Downloads/Normal Samples (tumor)/TCGA-B2-5641.tsv',
         '/Users/robertoruizfelix/Downloads/Tumor Samples/TCGA-CJ-5684.tsv',
         '/Users/robertoruizfelix/Downloads/Tumor Samples/TCGA-BP-4999.tsv',
         '/Users/robertoruizfelix/Downloads/Tumor Samples/TCGA-BP-4756.tsv']
for path in paths:
    # Extract the sample name from each path
    # [0] returns sample name as beginning of pathname
    sample_name = os.path.splitext(os.path.basename(path))[0]
    # Load new filename into pd df
    #sep (delimeter) by tabs bc tsv not csv
    #Skip first row as it is useless
    data = pd.read_csv(path, sep = '\t', skiprows = [0, 2, 3, 4, 5])    #skip unwanted rows
    data = data[['gene_id', 'fpkm_unstranded']]
    # Add sample name column
    data.columns = ['gene_id', sample_name]
    
    #Merge new data indo the existing df
    if merged_df.empty:
        merged_df = data
    else:
        # Includes all rows from both DataFrames and fills in 
#missing values with NaN where there are no matches
        # 'gene_id' is used to match rows
        merged_df = merged_df.merge(data, on = 'gene_id', how = 'outer')
print(data.columns)       # Show columns available

#Show first & last 5 rows
print(merged_df.head())

#Set 'gene_id' column as df anchor
#inplace = true bc we want to replace achor, not create a new df with a new anchor
merged_df.set_index('gene_id', inplace = True)
#Fill null values with 0's
merged_df.fillna(0, inplace = True)


#Calculting Fold Change

tumor_samples = ["TCGA-BP-4756", "TCGA-BP-4999", "TCGA-CJ-5684"]
normal_samples = ["TCGA-B2-5641", "TCGA-CW-5585", "TCGA-CZ-4863"]

#Compute mean expression value for each gene
tumor_mean = merged_df[tumor_samples].mean(axis = 1)
normal_mean = merged_df[normal_samples].mean(axis = 1)

# Compute Fold change 
fold_change = tumor_mean / normal_mean

# Replace +-inf with +-max, better data interpretation
max_fold_change = fold_change.replace([np.inf, -np.inf], np.nan).dropna().max()
fold_change.replace([np.inf], max_fold_change, inplace = True)
fold_change.replace([-np.inf], -max_fold_change, inplace = True)

#Extract upregulated (fc > 2) and downregulated (fc 0<x< 0.5) genes
upregulated = fold_change[fold_change > 2]
downregulated = fold_change[(0 < fold_change) & (fold_change < 0.5)]

#Create new df's for series of fold changes
# reset_index()resets the index of the series and converts it into a DataFrame
upregulated_df = upregulated.reset_index().rename(columns = {0: 'log2FoldChange', 'index': 'gene_id'})
downregulated_df = downregulated.reset_index().rename(columns = {0: 'log2FoldChange', 'index': 'gene_id'})

#Visualizing the data

#Exclude null values 
fold_change_graph = fold_change.dropna()
maroon = (128, 0, 0)
plt.hist(fold_change_graph, bins = 50, color = 'maroon', alpha = 0.7)
#A log2 scale, for example, allows you to easily see differences in expression that represent twofold changes. 
plt.yscale('log')
plt.xlabel('Fold change')
plt.ylabel('Number of Genes')
plt.title('Distribution of Fold Changes across Tumor Samples')
plt.show()

# Descriptive statistics 

fold_change_graph.describe()
print(f'Number of upregulated genes: {len(upregulated)}')
print(f'Number of downregulated genes: {len(downregulated)}')

#Save results (unique to each user)
upregulated_df.to_csv('Upregulated_genes_w_foldchange.csv', index = False)
downregulated_df.to_csv('Downregulated_genes_w_foldchange.csv', index = False)
output_path = '/Users/robertoruizfelix/Downloads/merged_Data.tsv'
merged_df.to_csv(output_path, sep = '\t', index = True)
print('file saved to: ', output_path)

