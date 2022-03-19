import numpy as np
import pandas as pd

from rpy2 import robjects
from rpy2.robjects import Formula

from rpy2.robjects import pandas2ri
pandas2ri.activate()

from rpy2.robjects.packages import importr

base = importr('base')
stats = importr('stats')
DESeq2 = importr('DESeq2')

df_genes = pd.read_csv('gene_counts.tsv', sep = '\t')
df_genes = df_genes.set_index('Gene')
df_genes = df_genes.T
df_genes.reset_index()

# just two colls
df_annot = pd.read_csv('annotation.csv')
df_annot = df_annot.loc[:,['Unnamed: 0','Biopsy Part']]
df_annot = df_annot.set_index('Unnamed: 0')

# concatination
df_res = pd.concat([df_genes, df_annot], axis = 1, join = 'inner')

#================================
#print(df_res['Biopsy Part'].unique()) # вывел ['Right' 'Left' 'Other' 'Rectum']
#print(df_annot['Biopsy Part'].unique()) # вывел ['Right' 'Left' 'Other' 'Rectum']
# У обоих dataframe"ов 562 строчки, пересечение верное
#================================

#Rectum сделам Left, удалим все образцы с Other
df_res.loc[df_res['Biopsy Part'] == 'Rectum','Biopsy Part'] = 'Left'
df_res = df_res.loc[df_res['Biopsy Part'].isin(['Right', 'Left'])] # убираю значение 'Other' тк в задании
                                                                   # сравнить экспрессию генов только между Left и Right

df_res_sorted = df_res.sort_values('Biopsy Part') # конечный dataframe
#print(df_res_sorted['Biopsy Part'].value_counts()) #Left = 269, Right = 176
df_res_sorted_final = df_res_sorted.iloc[:,:-1].T
print(df_res_sorted_final)

# здесь по сути забил на предыдущий датафрейм df_annot со столбцом Biopsy Part и сделал новый датафрейм
df_res_annot = pd.DataFrame({'BiopsyPart':['Left']*269 + ['Right']*176}, index=df_res_sorted_final.columns)
df_res_annot['BiopsyPart'] = stats.relevel(robjects.vectors.FactorVector(df_res_annot['BiopsyPart']), ref = 'Left') #это так и не понял что делает, взял из семинара
print(df_res_annot)

dds = DESeq2.DESeqDataSetFromMatrix(countData = df_res_sorted_final,
                                    colData = df_res_annot,
                                    design = Formula('~ BiopsyPart'))

dds = DESeq2.DESeq(dds)
res = DESeq2.results(dds, name = "BiopsyPart_Right_vs_Left")
res = DESeq2.lfcShrink(dds, coef = "BiopsyPart_Right_vs_Left", type = 'apeglm')
res = pd.DataFrame(base.as_data_frame(res))
res.index = df_res_sorted_final.index
res = res.sort_values('padj')
res = res.loc[res['padj'] < 0.05]
res = res.loc[res['log2FoldChange'].abs() >= 1] # вот здесь вопрос с тем, какой порог

res.to_csv('DESeq2_unpaired_left_vs_right_tummor_tissue.tsv', sep = '\t')
