import os                           # Get access to operating system
import pandas as pd                 # For data frame
import numpy as np                  # Calculations
import matplotlib.pyplot as plt     # Graph visualization

# Empty data frame to store data
merged_df = pd.DataFrame()

# List of working files, using pathnames
paths = []
for path in paths():
    # Extract the sample name from each path
    # [0] returns sample name as beginning of pathname
    filename = os.path.splitext(os.path.basename(path))[0]
    sample_name = '-'.join(filename.split('-'))[:3]
    # Load new filename into pd df
    #sep (delimeter) by tabs bc tsv not csv
    #Skip first row as it is useless
    data = pd.read_csv(path, sep = '\t', skiprows = 1)
    print(data.columns())       # Show columns available
for path in paths:      #Second optionary for loop just to show columns
    # Extract necessary columns
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

#Show first & last 5 rows
print(merged_df.head())
unwanted_rows = ['N_numapped', 'N_nultimapping', 'N_noFeature', 'N_ambiguous']
# ~ is the same as drop null values 
# 2 different ways to eliminate rows below
#merged_df = merged_df[merged_df['gene_id'].isin(unwanted_rows) == False]
merged_df = merged_df[~merged_df['gene_id'].isin(unwanted_rows)] 
#Set 'gene_id' column as df anchor
#inplace = true bc we want to replace achor, not create a new df with a new anchor
merged_df.set_index('gene_id', inplace = True)
#Fill null values with 0's
merged_df.fillna(0, inplace = True)


#Calculting Fold Change

tumor_samples = []
normal_samples = []

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
downregulated = fold_change[0 < fold_change < 0.5]

#Create new df's for series of fold changes
# reset_index()resets the index of the series and converts it into a DataFrame
upregulated_df = upregulated.reset_index().rename(columns = {0: 'log2DoldChange', 'index': 'gene_id'})
downregulated_df = downregulated.reset_index().rename(columns = {0: 'log2DoldChange', 'index': 'gene_id'})

#Visualizing the data

#Exclude null values 
fold_change_graph = fold_change.dropna()
maroon = (128, 0, 0)
plt.hist(fold_change_graph, bins = 50, color = maroon)
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
print(f'Top 5 upregulated Genes: ')
print(upregulated_df.sort_values(by = 'log2FoldChange', ascending = False).head())
print(f'\nTop 5 downregulated Genes: ')
print(downregulated_df.sort_values(by = 'log2FoldChange').head())

#Save results (unique to each user)
upregulated_df.to_csv('Upregulated_genes_w_foldchange.csv', index = False)
downregulated_df.to_csv('Downregulated_genes_w_foldchange.csv', index = False)
output_path = 'C: /Users/'

