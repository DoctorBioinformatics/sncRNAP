#!/usr/bin/env python
# coding: utf-8

# In[412]:


import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import venn
import warnings
from matplotlib.backends.backend_pdf import PdfPages
warnings.filterwarnings("ignore")


# ### Read ncRNA db 

# In[413]:


# open other ncRNA db file
fasta = sys.argv[1]
with open(fasta, "r") as file1:
    content1 = file1.read()

# open concatenated file and count ncRNA classes
lines = content1.splitlines()


# ### A) ncRNA length

# In[415]:


# iterate over ids and extract ncRNA classes
lengths = {}
ncRNA_id = None
for line in lines:
    if 'Mt' in line:
        continue

    if line.startswith(">"):
        if 'trna' in line:
            ncRNA_id = line.split('>')[1].strip() + '_tRNA'
        if 'MIMA' in line:
            ncRNA_id = line.split('>')[1].strip()
            ncRNA_id = line.split('_')[1].strip() + '_miRNA'
        if ('trna' not in line) and ('MIMAT' not in line):
            ncRNA_id = line.split('>')[1].strip()
        lengths[ncRNA_id] = 0
        
    if ('trna' not in line) and ('MIMAT' not in line) and ('ENSMUSG' not in line) and (ncRNA_id is not None):
        lengths[ncRNA_id] += len(line.strip())
                
    if ('trna' not in line) and ('MIMAT' not in line) and ('ENSMUSG' not in line):
        ncRNA_id = None

length_df = pd.DataFrame(list(lengths.items()), columns=['ncRNA_id', 'Length'])
length_df = length_df.dropna()
length_df.loc[length_df['ncRNA_id'].str.contains('tRNA'), 'Length'] = length_df.loc[length_df['ncRNA_id'].str.contains('tRNA'), 'Length'] - 100


# In[416]:


# create an empty list to store the results
data = []

# loop through each row in the length_df
for i, row in length_df.iterrows():
    # split the ncRNA_id column to extract the class
    ncRNA_class = row['ncRNA_id'].split('_')[1]
    length = row['Length']
    
    # add the ncRNA class and length to the data list
    data.append((ncRNA_class, length))

# create a new DataFrame from the data list
data_new = pd.DataFrame(data, columns=['ncRNA_class', 'Length'])
data_new
# remove length values > 500
data_new = data_new[data_new['Length'] <= 500]

# create a range of values to use for the bins
bins = range(0, int(data_new['Length'].max())+5, 5)

# bin the length values
data_new['binned'] = pd.cut(data_new['Length'], bins=bins)

# group the data by ncRNA class and binned length values
grouped = data_new.groupby(['ncRNA_class', 'binned']).count().reset_index()

# create a pivot table from the grouped data
pivot_table = grouped.pivot(index='ncRNA_class', columns='binned', values='Length')

# plot the heatmap and append the figure object to the list
fig1 = plt.figure(figsize=(30,5))
sns.heatmap(pivot_table, cmap='viridis',vmin=0, vmax=100, linewidth =1, square =True)
plt.xlabel('ncRNA Length (binned)', fontsize=12) 
plt.ylabel('ncRNA class', fontsize=12) 
plt.xticks(fontsize=10) 
plt.yticks(fontsize=10) 
plt.tight_layout() 

# ### B) ncRNA count

# In[417]:


# create a dictionary to store the counts ncRNA lasses
counts = {}

for line in lines:
    # ignore Mt_RNA
    if 'Mt' in line:
        continue
    
    # extract ncRNA class and ignore transcript ids
    if line.startswith(">"):
        if 'trna' in line:
            ncRNA_class = 'tRNA'
            counts[ncRNA_class] = counts.get(ncRNA_class, 0) + 1
            
        if 'MIMA' in line:
            ncRNA_class = 'miRNA'
            counts[ncRNA_class] = counts.get(ncRNA_class, 0) + 1
                   
        if ('trna' not in line) and ('MIMAT' not in line):
            ncRNA_class = line.split('>')[1].strip()
            ncRNA_class = line.split('_')[1].strip()
            counts[ncRNA_class] = counts.get(ncRNA_class, 0) + 1

# create a dataframe to store the result
result_df = pd.DataFrame(list(counts.items()), columns=['ncRNA_class', 'Count'])


# In[418]:


# plot gene count

# calculate the percentage of each class
result_df['Percentage'] = result_df['Count'] / result_df['Count'].sum() * 100

# Combine all classes with a percentage less than 5% into 'Others'
result_df.loc[result_df['Percentage'] < 5, 'ncRNA_class'] = 'Others'

# group the dataframe by the 'ncRNA_class' column
grouped_df = result_df.groupby('ncRNA_class').sum()

# plot the pie chart and append the figure object to the list
fig2, ax = plt.subplots(figsize=(8, 4))
ax.pie(grouped_df['Percentage'], labels=grouped_df.index, autopct='%1.1f%%',
       shadow=True, startangle=140)
ax.axis('equal')

# plot the bar chart and append the figure object to the list
fig3, ax2 = plt.subplots(figsize=(8, 4)) 
bar_counts = grouped_df['Count']
ax2.bar(grouped_df.index, grouped_df['Count'], alpha=0.3, color='gray')
ax2.set_ylim(0, grouped_df['Count'].max())
ax2.set_ylabel('Number of genes')
ax2.set_ylim(0, grouped_df['Count'].max())
for i, v in enumerate(bar_counts):
    ax2.text(i, v + 0.5, str(v), ha='center')


# ### C) ncRNA counts per sample

# In[419]:


# Read CSV file
counts_csv = sys.argv[2]
df = pd.read_csv(counts_csv)
df = df.drop(df.columns[0], axis=1)

# rename column id
df = df.rename(columns={df.columns[0]: 'id'})


# #### change id names to id_class

# In[420]:


# Iterate over the ids 
for i, row in enumerate(df['id']):
    
    if "trna" in row:    
        # change row name to row name + ncRNA_class
        ncRNA_class = 'tRNA'

        new_id = row + "_" + ncRNA_class
        df['id'] = df['id'].replace({row: new_id})
        
    if "MIMA" in row:
        ncRNA_class = 'miRNA'
        new_string = row.split('_')[1]
        new_id = new_string + "_" + ncRNA_class
        df['id'] = df['id'].replace({row: new_id})
    
    if ('trna' not in row) and ('MIMAT' not in row):
        new_id = row
        df['id'] = df['id'].replace({row: new_id})
        
# remove rows that don't have a classification
df = df[~df['id'].str.contains('Mt')]


# #### check ncRNA counts across samples

# In[421]:


# iterate over column samples and find ids that with values !=0 
counts = {}
for column_name, column_values in df.iloc[:, 4:].iteritems():
    sample_name = column_name
    sample_counts = {}
    for index, value in enumerate(column_values.values):
        if value == 0:
            continue
        row = df.iloc[index, 0]
        
        split_line = row.split('_')
        ncRNA_class = split_line[1].strip()
        sample_counts[ncRNA_class] = sample_counts.get(ncRNA_class, 0) + 1

        counts[sample_name] = sample_counts
    
# create a dataframe to store the result
genes_classes_df = pd.DataFrame(counts)

# barplot
sns.set_style("whitegrid")
sns.set_palette("colorblind")

fig4, ax = plt.subplots(figsize=(18, 6))

genes_classes_df.plot(kind='bar', stacked=False, ax=ax)

for i in range(genes_classes_df.shape[0]):
    for j in range(genes_classes_df.shape[1]):
        count = genes_classes_df.iloc[i, j]
        x_pos = i - 0.5 * (len(genes_classes_df.columns) - 1) / len(genes_classes_df.columns) + j / len(genes_classes_df.columns)
        y_pos = count + 5
        ax.text(x_pos, y_pos, str(count), ha='center', va='bottom')

ax.set_xlabel("ncRNA Classes")
ax.set_ylabel("Counts")
ax.set_title("ncRNA Counts for Each Sample")

ax.legend(title='Samples', frameon=False, loc='center left', bbox_to_anchor=(1.0, 0.5))


# #### venn diagram per ncRNA 

# In[422]:


# create a venn df and convert values to 0/1
venn_df = df.drop(df.columns[1:4], axis=1)

# convert values 
threshold = 0
venn_df.iloc[:, 1:] = venn_df.iloc[:, 1:].astype(float)
venn_df.iloc[:, 1:] = venn_df.iloc[:, 1:].applymap(lambda x: np.where(x > threshold, 1, 0))

# venn plot

# Extract sncRNA class from "id" column
venn_df["sncRNA_class"] = venn_df["id"].apply(lambda x: x.split("_")[1])

# Get unique classes
classes = venn_df["sncRNA_class"].unique()

# Get unique sample names
sample_names = [col for col in venn_df.columns if col not in ["id", "sncRNA_class"]]
sample_counts = len(sample_names)

# check if number of samples is less than or equal to 4
if sample_counts <= 4:
    # iterate through each class and plot Venn diagrams for each class
    venn_figs = [] # define a list to store all venn diagrams
    for class1 in classes:
        class1_data = venn_df[venn_df["sncRNA_class"] == class1].iloc[:, 1:]
        sets = []
        for sample in sample_names:
            sets.append(set(class1_data.index[class1_data[sample] > 0]))

        venn_function = getattr(venn, f'venn{sample_counts}')
        labels = venn.get_labels(sets, fill=['number'])
        fig5, ax = venn_function(labels, names=sample_names)
        ax.set_title(f"{class1} in {', '.join(sample_names)}")

        venn_figs.append(fig5)
else:
    venn_figs = []
    print("Number of samples is greater than 4. Venn plot cannot be generated.")



# #### expression profile per sample

# In[423]:


# Create bins for expression values
bins = [0, 1, 10, 100, 1000, float("inf")]
labels = ["<=1", "1-10", "10-100", "100-1000", ">1000"]

# Create a list to store the binned values for each sample
binned_values = []
for column_name, column_values in df.iloc[:, 4:].iteritems():
    sample_name = column_name
    
    # Bin the expression values
    binned_expression_values = pd.cut(column_values, bins=bins, labels=labels, right=False, include_lowest=True)
    
    # Store the binned expression values for each sample
    binned_values.append(pd.DataFrame({'Expression Level': binned_expression_values, 'Sample': sample_name}))

# Concatenate the binned values for each sample
binned_df = pd.concat(binned_values)

# Pivot the data so that each expression level is a column and each sample is a row
binned_df = binned_df.pivot_table(index='Sample', columns='Expression Level', aggfunc='size').fillna(0)

# Plot the binned values as a stacked bar plot
# fig6, ax = plt.subplots()
fig6, ax = plt.subplots(figsize=(18, 6))

binned_df.plot(kind='barh', stacked= True, ax=ax)
ax.set_xlabel("Number of ncRNAs")
ax.set_ylabel("Sample")
ax.set_title("ncRNA Expression Bins")
ax.legend(title="Expression (RPM)",loc="center left", bbox_to_anchor=(1.0, 0.5))
# plt.show()


# #### expression profile per ncRNA class

# In[424]:


# Create bins for expression values
bins = [0, 1, 10, 100, 1000, float("inf")]
labels = ["<=1", "1-10", "10-100", "100-1000", ">1000"]

# Create a list to store the binned values for each sample and ncRNA class
binned_values = []
for column_name, column_values in df.iloc[:, 4:].iteritems():
    sample_name = column_name
    
    # Bin the expression values
    binned_expression_values = pd.cut(column_values, bins=bins, labels=labels, right=False, include_lowest=True)
    
    # Store the binned expression values for each sample and ncRNA class
    sample_df = pd.DataFrame({'Expression Level': binned_expression_values, 'Sample': sample_name})
    sample_df['ncRNA Class'] = df['id'].apply(lambda x: x.split('_')[1])
    binned_values.append(sample_df)

# Concatenate the binned values for each sample and ncRNA class
binned_df = pd.concat(binned_values)

# Plot each ncRNA class in a separate facet plot
g = sns.FacetGrid(binned_df, col='ncRNA Class', col_wrap=3, height=4, aspect=2, sharex=False, sharey=False, gridspec_kws={'wspace': 1, 'hspace': 1})
g.map_dataframe(sns.histplot, x='Sample', hue='Expression Level', multiple='stack')
g.set_axis_labels('Sample', 'Number of ncRNAs')
g.add_legend(title='Expression Level', bbox_to_anchor=(1, 0.5), loc='right', borderaxespad=0., ncol=1, title_fontsize=14, labels=labels[::-1])
g.set_xticklabels(rotation=90, fontsize=6)
fig7 = g.fig


# ### D) DE ncRNAs

# In[425]:


df.info()


# #### DE gene counts

# In[426]:


# plot log2fc distribution
fig8 = plt.figure()
plt.hist(df['log2FoldChange'], bins=50, edgecolor='black')
plt.xlabel('log2FoldChange')
plt.ylabel('Frequency')
plt.title('Histogram of log2FoldChange values for all ncRNAs')

# plot padj distribution
fig9 = plt.figure()
plt.hist(df['padj'], bins=50, edgecolor='black')
plt.xlabel('padj')
plt.ylabel('Frequency')
plt.title('Histogram of padj values for all ncRNAs')

# determine DE value threshold nad plot percentage 
counts = ((df['log2FoldChange'] > 1) | (df['log2FoldChange'] < -1)) & (df['padj'] < 0.005)
count = counts.sum()

# Plot the pie chart
fig10 = plt.figure()
labels = ['Differentially Expressed Genes', 'Non-Differentially Expressed Genes']
sizes = [count, len(df)]
ax1 = fig10.add_subplot(111)
ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
        shadow=True, startangle=90)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle
ax1.set_title(str(count)+' DE genes')


# #### DE gene counts per class

# In[427]:


# count DE ncRNAs
differentially_expressed = df[((df['log2FoldChange'] < -1) | (df['log2FoldChange'] > 1)) & (df['padj'] < 0.005)]

# create a dictionary to store the counts ncRNA lasses
counts = {}

# iterate over the rows in the df2 gene id
for row in differentially_expressed[differentially_expressed.columns[0]]:
    split_line = row.split('_')
    ncRNA_class = split_line[1].strip()
    counts[ncRNA_class] = counts.get(ncRNA_class, 0) + 1

# create a dataframe to store the result
DE_classes_df = pd.DataFrame(list(counts.items()), columns=['ncRNA_class', 'Count'])
DE_classes_df

# Get the count data
counts = DE_classes_df['Count'].tolist()
labels = DE_classes_df['ncRNA_class'].tolist()

# Plot the bar plot using seaborn
fig11 = plt.figure()
sns.barplot(x='ncRNA_class', y='Count', data=DE_classes_df)

# Add the count on top of each bar
for i in range(len(DE_classes_df)):
    count = DE_classes_df['Count'].iloc[i]
    plt.text(i, count, count, ha='center', va='bottom')

plt.ylabel("Number of DE genes")
# plt.show()


# Plot the pie chart
fig12 = plt.figure()
plt.pie(counts, labels=labels, startangle=90, counterclock=False, autopct='%1.1f%%')
plt.axis('equal')
plt.title("Number of DE genes")
# plt.show()


# ####  padj, log2fc, and DE gene % for each class

# In[428]:


# Split the id column to extract ncRNA class
df['ncRNA_class'] = df['id'].apply(lambda x: x.split("_")[1])

# Define a function to plot the histogram of log2FoldChange and pvalue for each ncRNA class
def plot_histogram(feature, grid):
    grid.map(sns.histplot, feature, bins=30)

# Define a function to plot a pie chart for each ncRNA class
def plot_piechart(data, **kwargs):
    counts = ((data['log2FoldChange'] > 1) | (data['log2FoldChange'] < -1)) & (data['padj'] < 0.005)
    count = counts.sum()
    labels = ['DE', 'Non-DE']
    sizes = [count, len(data) - count]
        
    plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
    plt.axis('equal')

# Plot histograms for log2fc and pvalues
g1 = sns.FacetGrid(df, col='ncRNA_class', col_wrap=3, sharex=False, sharey=False)
fig13 = g1.fig
plot_histogram('log2FoldChange', g1)

g2 = sns.FacetGrid(df, col='ncRNA_class', col_wrap=3, sharex=False, sharey=False)
fig14 = g2.fig
plot_histogram('padj', g2)

# Plot a pie chart for each ncRNA class
g3 = sns.FacetGrid(df, col='ncRNA_class', col_wrap=3, sharex=False, sharey=False, height=3)
g3.map_dataframe(plot_piechart)
fig15 = g3.fig


# #### volcano plot per class

# In[429]:


# Volcano plot
def plot_volcano(filtered_df, de_genes_class, ncRNA_class, ax=None):
    if ax is None:
        ax = plt.gca()
    # Plot the differentially expressed genes in red
    ax.scatter(de_genes_class["log2FoldChange"], -np.log10(de_genes_class["pvalue"]), color='red', label='DE genes')

    # Add the id labels for the DE genes
    for i, row in de_genes_class.iterrows():
        ax.annotate(row["id"], (row["log2FoldChange"], -np.log10(row["pvalue"])), fontsize=7)

    # Plot the remaining genes in blue
    non_de_genes_class = filtered_df[~filtered_df.index.isin(de_genes_class.index)]
    ax.scatter(non_de_genes_class["log2FoldChange"], -np.log10(non_de_genes_class["pvalue"]), color='blue', label='non-DE genes')

    # Add the ncRNA label to the plot
    ax.set_xlabel("log2FoldChange")
    ax.set_ylabel("-log10(pvalue)")
    ax.set_title(f"Volcano plot for {ncRNA_class} ncRNA class")
    ax.legend()

# Get the list of unique ncRNA classes from the "id" column
ncRNA_classes = df["id"].apply(lambda x: x.split("_")[1]).unique()

volcano_list = []
for ncRNA_class in ncRNA_classes:
    # Filter the data for the current ncRNA class
    filtered_df = df[df["id"].apply(lambda x: x.split("_")[1] == ncRNA_class)]
    de_genes_class = differentially_expressed[differentially_expressed["id"].apply(lambda x: x.split("_")[1] == ncRNA_class)]
    fig16, ax = plt.subplots(figsize=(10, 10))
    plot_volcano(filtered_df, de_genes_class, ncRNA_class, ax)
    volcano_list.append(fig16)
    


# #### heatmap for DE genes

# In[430]:


# Subset the differentially expressed genes
diff_exp_gene_values = differentially_expressed.iloc[:,4:]

# Create the heatmap figure
fig17, ax = plt.subplots(figsize=(22,14))
sns.heatmap(diff_exp_gene_values, xticklabels=diff_exp_gene_values.columns, yticklabels=differentially_expressed["id"], cmap="YlGnBu", annot=True, annot_kws={"size": 10}, ax=ax)

# Add labels and title
plt.xlabel("Sample", fontsize=12)
plt.ylabel("ID", fontsize=12)
plt.title("Expression (RPM)", fontsize=15)
plt.xticks(rotation=45, ha="right", fontsize=10)
plt.yticks(fontsize=10)

# Adjust the margins to show all the IDs
plt.subplots_adjust(left=0.25, bottom=0.25)

# # Show the plot
# plt.show()


# ### Export figs as combined PDF

# In[431]:


# create a list of the figures to be saved
figures = [fig1, fig2, fig3, fig4] + [fig6, fig7, fig8, fig9, fig10, fig11, fig12, fig13, fig14, fig15] + volcano_list + [fig17]

if venn_figs:
    figures += venn_figs

# create a PDF file using PdfPages
with PdfPages('sncRNA_profiling.pdf') as pdf:
    for fig in figures:
        pdf.savefig(fig)


# ### Export log2fc values and gene sequences

# #### log2fc values for top candidates

# In[432]:


differentially_expressed


# In[433]:


# save to csv 
differentially_expressed.to_csv('de_genes.csv', index=False)


# In[434]:


differentially_expressed


# #### ncRNA sequence for top candidates

# In[445]:


# extract the 'id' column as a list
ids = differentially_expressed['id'].str.split('_').str[0].tolist()

# initialize an empty dictionary to store the sequences
sequences = {}

# initialize an empty string to store the current sequence
sequence = ''

# initialize a flag to keep track of header lines
header_line = True

for id in ids:
    
    # iterate over the lines in the fasta file
    for line in lines:
        # if the line starts with a '>', it is a header line
        if line.startswith('>'):
            # if it's not the first header line
            if not header_line:
                # check if the header is in the 'ids' list
                if id in header:
                    # store the sequence in the dictionary
                    sequences[header] = sequence
            # remove the '>' symbol and split the header on spaces
            header = line[1:].strip().split()[0]
            # reset the sequence
            sequence = ''
            # update the header line flag
            header_line = True
        # if the line does not start with a '>', it is a sequence line
        else:
            # add the line to the current sequence
            sequence += line.strip()
            # update the header line flag
            header_line = False

# # store the final sequence in the dictionary
sequences[header] = sequence

# open the output file for writing
with open('top_candidates.fa', 'w') as outfile:
    for id in ids:
        for seq in sequences:
            if id in seq:
                # write the sequences for each id to the output file
                print(f">{seq}\n{sequences[seq]}\n")
                outfile.write(f">{seq}\n{sequences[seq]}\n")

print("Sequences saved to top_candidates.fa")
print("Analysis completed.. Check summary PDF, candidate fasta sequences, and de_genes.csv")
