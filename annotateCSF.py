#!/usr/bin/env python

#tk imports
from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
from tkinter import simpledialog
from tkinter.colorchooser import askcolor
import PIL.Image
import PIL.ImageTk
#external imports
import os
import sys
import pandas as pd
import scvi as scvi
import scanpy as sc
import skmisc
import matplotlib.pyplot as plt
import seaborn as sns
import scrublet as scr
import torch
import numpy as np
import random
import statsmodels.api as sm
from statsmodels.formula.api import ols
import pandastable as pt
from termcolor import colored
import termcolor

# root for main widget
root = Tk()
#root.iconbitmap("files/icon3.ico")
root.title('annotate CSF')

# welcome statement
def welcome_statement():
    print(colored("\nWelcome to aCSF! \nYou can use aCSF for label transfer and subsequent plotting and analysis of your CSF scRNAseq data. Please refer to the reference manual for more detailed information regarding the usage of aCSF.\n", 'yellow'), colored("\nThis tool was built using scanpy, anndata, scVI, and scANVI workflows. When using, please consider citing the original publications as well: \n\nWolf et al. (2018), Scanpy: large-scale single-cell gene expression data analysis, Genome Biology. \n\nVirshup et al. (2021) anndata: Annotated data, bioRxiv. \n\nGayoso, Lopez, Xing, et al (2022), A Python library for probabilistic analysis of single-cell omics data, Nature Biotechnology. \n\nLopez (2018), Deep generative modeling for single-cell transcriptomics, Nature Methods. \n\nXu, Lopez, Mehlman, et al (2021), Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models, Molecular Systems Biology.\n-----------------------------------------------------------------------------------------------------------------------\n", "white" ))
welcome_statement()

# get logo
if getattr(sys, 'frozen', False):
    image = PIL.Image.open(os.path.join(sys._MEIPASS, "files/acsf_logo.png"))
else:
    image = PIL.Image.open("files/acsf_logo.png")

image = image.resize((1450,560))
image =  PIL.ImageTk.PhotoImage(image)
wl = Label(root, image = image)

# get directory for 10x data
def get_path():
    path = filedialog.askdirectory(title = 'Select folder with 10X-formatted files')
    path_var.set(value = path)
    print(colored('Data will be loaded from:'))
    print(path_var.get())

# get directory to store .h5ad
def get_save_path():
    path = filedialog.askdirectory(title = 'Select folder to store your file')
    path_var.set(value = path)
    print(colored('Data will be saved to:'))
    print(path_var.get())

# get genes or metadata to plot in UMAP	
def get_genes_to_plot():
    genes = simpledialog.askstring(title='Genes to plot', prompt='Name gene or metadata to plot, example1: FOXP3; example2: condition')
    gene_var.set(value=genes)
    global genes_to_plot
    genes_to_plot = gene_var
# get genes to plot in heatmap
def get_genes_to_plot2():
    genes = simpledialog.askstring(title='Genes to plot', prompt='Please name the desired genes seperated by commas, e. g. "CCR6, FOXP3, CCR5"')
    gene_var.set(value=genes)
    global genes_to_plot
    genes_to_plot = gene_var

# get order of conditions for box-/jitter plot from user
def get_condition_order():
    print(colored('Available conditions in your dataset:', 'yellow'))
    print(bdata.obs.condition.value_counts().index.tolist())
    order = simpledialog.askstring(title='Change Condition order',  prompt='Want to change condition order in plot? Type new order as e.g.: \'Control\', \'Disease1\', \'Disease2\'')
    order_var.set(value=order)
    global condition_order
    condition_order = order_var

# get hues to color conditions accordingly in box-/jitter plot	
def get_hues():
    cols = {}
    col_prompts = list(['Choose first color', 'Choose second color', 'Choose third color', 'Choose fourth color', 'Choose fifth color', 'Choose sixth color', 'Choose seventh color', 'Choose eighth color', 'Choose nineth color', 'Choose tenth'])
    print(len(bdata.obs.condition.value_counts().index.tolist()))
    for i in range(0, len(bdata.obs.condition.value_counts().index.tolist())):
        col_choice = askcolor(title=col_prompts[i])
        cols[i] = str(col_choice[1])
    hue = list(cols.values())
    hue_var.set(value=hue)
    global new_hues
    new_hues = hue_var

	
def plot_heatmap_def():
    plot1 = {' B activated | B IL4R+ | B atypical | Plasmacells': ['CD19', 'CD79A', 'TNFRSF13B'], 'pDCs': ['IL3RA', 'IRF8', 'LAMP5'], 'Monocytes | BAM MRC1+ | BAM EMP3+ | MG CX3CR1+ | MG CCL2+ | MG TREM2hi': ['CD14', 'CD68', 'MS4A7'], 'mDCs CD1c+ | mDCs AXL+SIGLEC6+ | mDCs CLEC9A+': ['CD1C', 'AFF3', 'HLA-DRB5'], 'NK bright | NK dim | TR-NK | ILC': ['NCAM1', 'GNLY', 'XCL1'], 'MAIT | gdT Vd2+ | gdT Vd2- | CD8 CM | CD8 EM HLA-DRA+ | CD8 EM CD160+ | CD8 TRM ITGA1+ | CD8 TRM ITGA1- | CD8 CTL': ['CD8B', 'CD8A', 'GZMH'], 'Tfh | Th17 | Th2/Th22 | Th1 | Tregs | CCR5high Th17.1 | CD4 TEMRA': ['CD4', 'TNFRSF25', 'AQP3']}
    plot2 = {' B activated': ['CD69', 'CXCR4', 'BACH2'], 'B IL4R+': ['IL4R', 'FCER2', 'IGHM'], 'B atypical': ['FCRL5', 'CD1C', 'GPR34'], 'Plasmacells': ['IGHG1', 'SDC1', 'CD38'], 'pDCs': ['IL3RA', 'IRF7', 'IRF8'], 'Monocytes': ['VCAN', 'S100A12', 'CCR2'], 'BAM MRC1+': ['MRC1', 'KCNAB1', 'MARCO'], 'BAM EMP3+': ['EMP3', 'CYP27A1', 'TIMD4'], 'MG CX3CR1+': ['CX3CR1', 'TMEM119', 'P2RY12'], 'MG CCL2+': ['SPP1', 'IRAK2', 'ITGAX'], 'MG TREM2hi': ['TREM2', 'APOC1', 'GPNMB'], 'mDCs CD1c+': ['CD1C', 'CD207', 'CD1E'], 'mDCs AXL+SIGLEC6+': ['SIGLEC6', 'CLIC3', 'CD5D'], 'mDCs CLEC9A+': ['CLEC9A', 'PPP1R14A', 'TMEM14A'], 'NK bright': ['NCAM1', 'KLRC3', 'SPINK2'], 'NK dim': ['EOMES', 'TTC38', 'FGFBP2'], 'TR-NK': ['IKZF3', 'KRT81', 'ITGAE'], 'ILC': ['KIT', 'TNFRSF4', 'IFNGR2'], 'MAIT': ['TRAV1-2', 'NCR3', 'KLRB1'], 'gdT Vd2+': ['TRDC', 'TRGV9', 'ZBTB16'], 'gdT Vd2-': ['TRDV2', 'LSR', 'RTKN2'], 'CD8 CM': ['CCR7', 'NELL2', 'MAL'], 'CD8 EM HLA-DRA+': ['HLA-DRA', 'MSC', 'HLA-DRB5'], 'CD8 EM CD160+': ['CD160', 'FCRL6', 'FGR'], 'CD8 TRM ITGA1+': ['ITGA1', 'ZNF683', 'ITGAE'], 'CD8 TRM ITGA1-': ['NR4A2', 'CEMIP2', 'CXCR4'], 'CD8 CTL': ['GNLY', 'FXYD2', 'GZMB'], 'Th17.1,CD4': ['TEMRA', 'PASK', 'CXCR5'], 'Th17': ['CCR6', 'USP10', 'RORC'], 'Th2/Th22': ['SLC40A1', 'SOCS1', 'CCR4'], 'Th1': ['TBX21', 'TRBC1', 'CXCR3'], 'Tregs': ['FOXP3', 'RTKN2', 'IKZF2'], 'CCR5high Th17.1': ['CXCR6', 'CD69', 'RGS1'], 'CD4 TEMRA': ['GZMH', 'PDCD1', 'CX3CR1']}

    subsets_in_data = set(adata_sub2.obs.predictions.value_counts().index.values)   # gibt eine Liste der vorhandenen subsets im dataset
    genes_in_data_set = set(adata_sub2.var.gene_symbols)   # extrahiert die Liste von vorhandenen Genen aus dem dataset

    nr_subsets_in_data = 0
    for k in plot1.keys():
        if (len(subsets_in_data.intersection(k.split(' | '))) > 0):
            nr_subsets_in_data += 1

    # plot lineage heatmap
    if (nr_subsets_in_data > 1):
        print ('Generating lineage heatmap...')
        genes_in_data = genes_in_data_set.intersection([item for sublist in plot1.values() for item in sublist])
        plot_data = {}
        for k,v in plot1.items():
            gene = next((gene for gene in v if (gene in genes_in_data)), None)
            for subset in k.split(' | '):
                if (subset in subsets_in_data and gene is not None):
                    plot_data[subset] = gene
                else:
                    print ("Omitting subset %s due to missingness of the subset or lack of marker gene expression." % subset)
        sc.pl.matrixplot(adata_sub2[adata_sub2.obs.predictions.isin(list(plot_data.keys()))], categories_order=list(plot_data.keys()), var_names=list(dict.fromkeys(plot_data.values())), groupby='predictions', standard_scale='var', cmap='inferno')

    # plot subset heatmap
    print ('Generating subset heatmap...')
    genes_in_data = genes_in_data_set.intersection([item for sublist in plot2.values() for item in sublist])
    plot_data = {}
    for k,v in plot2.items():
        gene = next((gene for gene in v if (gene in genes_in_data)), None)
        if (gene is not None and k in subsets_in_data):
            plot_data[k] = gene
        else:
            print ("Omitting subset %s due to missingness of the subset or lack of marker gene expression." % k)
    sc.pl.matrixplot(adata_sub2[adata_sub2.obs.predictions.isin(list(plot_data.keys()))], categories_order=list(plot_data.keys()), var_names=list(plot_data.values()), groupby='predictions', standard_scale='var', cmap='inferno')

def plot_heatmap_cust():
    global genes_to_plot
    get_genes_to_plot2()
    genes_to_plot = genes_to_plot.get()

    print ('Generating custom heatmap...')
    genes_in_data_set = set(adata_sub2.var.gene_symbols)   # extrahiert die Liste von vorhandenen Genen aus dem dataset
    plot_genes = genes_in_data_set.intersection(genes_to_plot.split(', '))

    sc.pl.matrixplot(adata_sub2, var_names=list(plot_genes), groupby='predictions', standard_scale='var', dendrogram=True, cmap='inferno')

def plot_heatmap():
    if (heatmap_choice.get() == 'Custom heatmap'):
        plot_heatmap_cust()
    else:
        plot_heatmap_def()

def plot_umap():
    subset = subset_choice2.get()
    global genes_to_plot
    genes_to_plot = genes_to_plot.get()

    if ((genes_to_plot not in ['predictions', 'orig.ident', 'condition', 'study']) and (genes_to_plot not in adata_to_sub.var.gene_symbols)):
        print('You entered ' + genes_to_plot + '. This is either not a valid gene symbol or the gene was not found in the dataset')
        return

    if subset == 'Do not subset':
        try:
            sc.pl.umap(adata_to_sub, color = [genes_to_plot], size=40, cmap="inferno", legend_loc="on data")
        except:
            print('You entered ' + genes_to_plot + '. This is either not a valid gene symbol or the gene was not found in the dataset')

    print("subset chosen is:" + subset)
    if subset == 'CD4 T cells':
        sub_to_keep = ['Tfh', 'Th17', 'Th2/Th22', 'Th1', 'Tregs', 'CCR5high Th17.1', 'CD4 TEMRA']
    elif subset == 'CD8 T cells':
        sub_to_keep = ['CD8 TRM ITGA1+', 'CD8 TRM ITGA1-', 'CD8 EM HLA-DRA+', 'CD8 EM CD160+', 'CD8 CM', 'CD8 CTL', 'gdT Vd2+', 'gdT Vd2-', 'MAIT']
    elif subset == 'Myeloid cells':
        sub_to_keep = ['MG CX3CR1+', 'Monocytes', 'MG CCL2+', 'MG TREM2hi', 'BAM MRC1+', 'BAM EMP3+', 'mDCs CD1c+', 'mDCs AXL+SIGLEC6+', 'mDCs CLEC9A+']
    elif subset == 'NK cells':
        sub_to_keep = ['NK bright', 'NK dim', 'TR-NK', 'ILC']
    elif subset == 'B cells':
        sub_to_keep = ['B IL4R+', ' B activated', 'Plasmacells', 'B atypical']

    if (subset != 'Do not subset') and ('adata_sub' not in globals()):
        print(colored('You chose to subset your data for the UMAP plot. To provide better resolution, a basic scanpy workflow will be run for your chosen subset.', 'yellow'))
        global adata_sub
        adata_sub = bdata.copy()
        if 'Th2/Th22' in adata_sub.obs.predictions.values:
             adata_sub.obs['mnc_lineage'] = 'CD4 T cells'
        elif 'CD8 EM HLA-DRA+' in adata_sub.obs.predictions.values:
             adata_sub.obs['mnc_lineage'] = 'CD8 T cells'
        elif 'BAM EMP3+' in adata_sub.obs.predictions.values:
            adata_sub.obs['mnc_lineage'] = 'Myeloid cells'
        elif 'NK bright' in adata_sub.obs.predictions.values:
            adata_sub.obs['mnc_lineage'] = 'NK cells'
        elif ' B activated' in adata_sub.obs.predictions.values:
            adata_sub.obs['mnc_lineage'] = 'B cells'
        print(adata_sub.obs['mnc_lineage'].value_counts())
        adata_sub = adata_sub[adata_sub.obs.predictions.isin(sub_to_keep)].copy()
        sc.pp.highly_variable_genes(adata_sub, n_top_genes=1500, subset=True, flavor = 'seurat', batch_key='orig.ident')
        answer = messagebox.askyesno(title="scanpy or scVI integration?", message = "Click yes for scVI workflow (slower, but probably better integration), click 'No' for scanpy workflow (faster)")
        if answer == True:
            scvi.data.setup_anndata(adata_sub, layer="counts", batch_key="orig.ident", categorical_covariate_keys = ['tissue', 'study'], continuous_covariate_keys = ['pct_counts_mt'], labels_key='labels')
            print(colored('Running scVI workflow. Commencing training...', 'red'))
            sub_model = scvi.model.SCVI(adata_sub, n_hidden = 128, n_latent = 10)
            sub_model.train(max_epochs = 400, early_stopping=True, early_stopping_patience = 15)
            print(colored('...done!', 'cyan'))
            adata_sub.obsm['X_scVI'] = sub_model.get_latent_representation()
            sc.pp.neighbors(adata_sub, use_rep = 'X_scVI')
        else:
            print(colored('Running scanpy workflow for subset. Batch effect correction via BBKNN...', 'red'))
            sc.pp.regress_out(adata_sub, ['total_counts', 'pct_counts_mt'])
            sc.pp.scale(adata_sub)
            sc.tl.pca(adata_sub)
            sc.pp.neighbors(adata_sub, n_pcs = 30, n_neighbors = 10)
            sc.external.pp.bbknn(adata_sub, batch_key='orig.ident')
        print(colored('...done!', 'cyan'))
        print(colored('Computing UMAP. Leiden clustering will also be applied at standard resolution (1.0), you can view this by typing: "leiden" in the dialogue IV. a', 'red'))
        sc.tl.umap(adata_sub)
        sc.tl.leiden(adata_sub, resolution = 1.0)
        print(colored('...done!', 'cyan'))
        try:
            sc.pl.umap(adata_sub, color=[genes_to_plot], size=40, legend_loc="on data")
        except:
            print('You entered ' + genes_to_plot + '. This is either not a valid gene symbol or the gene was not found in the dataset')
    else:
        if (subset != 'Do not subset') and ('adata_sub' in globals()) and (subset in adata_sub.obs.mnc_lineage.values):
            try:
                sc.pl.umap(adata_sub, color=[genes_to_plot], size=40, legend_loc="on data")
            except:
                print('You entered ' + genes_to_plot + '. This is either not a valid gene symbol or the gene was not found in the dataset')

    if 'adata_sub' in globals():
        if 'Th2/Th22' in adata_sub.obs.predictions.values:
             adata_sub.obs['mnc_lineage'] = 'CD4 T cells'
        elif 'CD8 EM HLA-DRA+' in adata_sub.obs.predictions.values:
             adata_sub.obs['mnc_lineage'] = 'CD8 T cells'
        elif 'BAM EMP3+' in adata_sub.obs.predictions.values:
            adata_sub.obs['mnc_lineage'] = 'Myeloid cells'
        elif 'NK bright' in adata_sub.obs.predictions.values:
            adata_sub.obs['mnc_lineage'] = 'NK cells'
        elif ' B activated' in adata_sub.obs.predictions.values:
            adata_sub.obs['mnc_lineage'] = 'B cells'

    global sub_to_keep2
    if (subset != 'Do not subset') and ('adata_sub' in globals() and (subset not in adata_sub.obs.mnc_lineage.values)):
        print(colored('You chose to subset your data for the UMAP plot. To provide better resolution, a model will be trained for your chosen subset.', 'yellow'))
        if subset == 'CD4 T cells':
            sub_to_keep2 = ['Tfh', 'Th17', 'Th2/Th22', 'Th1', 'Tregs', 'CCR5high Th17.1', 'CD4 TEMRA']
        elif subset == 'CD8 T cells':
            sub_to_keep2 = ['CD8 TRM ITGA1+', 'CD8 TRM ITGA1-', 'CD8 EM HLA-DRA+', 'CD8 EM CD160+', 'CD8 CM', 'CD8 CTL', 'gdT Vd2+', 'gdT Vd2-', 'MAIT']
        elif subset == 'Myeloid cells':
            sub_to_keep2 = ['MG CX3CR1+', 'Monocytes', 'MG CCL2+', 'MG TREM2hi', 'BAM MRC1+', 'BAM EMP3+', 'mDCs CD1c+', 'mDCs AXL+SIGLEC6+', 'mDCs CLEC9A+']
        elif subset == 'NK cells':
            sub_to_keep2 = ['NK bright', 'NK dim', 'TR-NK', 'ILC']
        elif subset == 'B cells':
            sub_to_keep2 = ['B IL4R+', ' B activated', 'Plasmacells', 'B atypical']
        mnc_lin = adata_sub.obs['mnc_lineage']
        adata_sub = bdata.copy()
        adata_sub.obs['mnc_lineage'] = mnc_lin
        adata_sub = adata_sub[adata_sub.obs.predictions.isin(sub_to_keep2)].copy()
        sc.pp.highly_variable_genes(adata_sub, n_top_genes=1500, subset=True, flavor = 'seurat', batch_key='orig.ident')
        answer = messagebox.askyesno(title="scanpy or scVI integration?", message = "Click yes for scVI workflow (slower, but probably better integration), click 'No' for scanpy workflow (faster)")
        if answer == True:
            scvi.data.setup_anndata(adata_sub, layer="counts", batch_key="orig.ident", categorical_covariate_keys = ['tissue', 'study', 'condition'], continuous_covariate_keys = ['pct_counts_mt'], labels_key='labels')
            print(colored('Running scVI workflow. Commencing training...', 'red'))
            sub_model = scvi.model.SCVI(adata_sub, n_hidden = 128, n_latent = 10)
            sub_model.train(max_epochs = 400, early_stopping=True, early_stopping_patience = 15)
            print(colored('...done!', 'cyan'))
            adata_sub.obsm['X_scVI'] = sub_model.get_latent_representation()
            sc.pp.neighbors(adata_sub, use_rep = 'X_scVI')
        else:
            print(colored('Running scanpy workflow for subset. Batch effect correction via BBKNN...', 'red'))
            sc.pp.regress_out(adata_sub, ['total_counts', 'pct_counts_mt'])
            sc.pp.scale(adata_sub)
            sc.tl.pca(adata_sub)
            sc.pp.neighbors(adata_sub, n_pcs = 30, n_neighbors = 10)
            sc.external.pp.bbknn(adata_sub, batch_key='orig.ident')
        print(colored('...done!', 'cyan'))
        print(colored('Computing UMAP. Leiden clustering will also be applied at standard resolution (1.0), you can view this by typing: "leiden" in the dialogue IV. a', 'red'))
        sc.tl.umap(adata_sub)
        sc.tl.leiden(adata_sub, resolution = 1.0)
        print(colored('...done!', 'cyan'))
        if ((genes_to_plot not in ['predictions', 'orig.ident', 'condition', 'study']) and (genes_to_plot not in adata_sub.var.gene_symbols)):
            sc.pl.umap(adata_sub, color=[genes_to_plot], size=40, legend_loc="on data")

# subset to mnc
def keep_mncs():
    global adata_sub2
    adata_sub2 = bdata.copy()

# subset function
def subset_fun():
    global adata_sub2
    adata_to_sub = bdata.copy()
    subset = subset_choice.get()
    if subset == 'CD4 T cells':
        adata_sub2 = adata_to_sub[adata_to_sub.obs.predictions.isin(['Tfh', 'Th17', 'Th2/Th22', 'Th1', 'Tregs', 'CCR5high Th17.1', 'CD4 TEMRA'])]
    elif subset == 'CD8 T cells':
        adata_sub2 = adata_to_sub[adata_to_sub.obs.predictions.isin(['CD8 TRM ITGA1+', 'CD8 TRM ITGA1-', 'CD8 EM HLA-DRA+', 'CD8 EM CD160+', 'CD8 CM', 'CD8 CTL', 'gdT Vd2+', 'gdT Vd2-', 'MAIT'])]
    elif subset == 'Myeloid cells':
        adata_sub2 = adata_to_sub[adata_to_sub.obs.predictions.isin(['MG CX3CR1+', 'Monocytes', 'MG CCL2+', 'MG TREM2hi', 'BAM MRC1+', 'BAM EMP3+', 'mDCs CD1c+', 'mDCs AXL+SIGLEC6+', 'mDCs CLEC9A+'])]
    elif subset == 'NK cells':
        adata_sub2 = adata_to_sub[adata_to_sub.obs.predictions.isin(['NK bright', 'NK dim', 'TR-NK', 'ILC'])]
    elif subset == 'B cells':
        adata_sub2 = adata_to_sub[adata_to_sub.obs.predictions.isin(['B IL4R+', ' B activated', 'Plasmacells', 'B atypical'])]
    elif subset == 'Do not subset':
        adata_sub2 = adata_to_sub.copy()

# recall full anndata after subsetting
def recall_full_data():
    global adata_to_sub
    adata_to_sub = bdata.copy()

def plot_freq1():
    anot_df = adata_sub2.obs['predictions'].value_counts()
    #anot_df = pbmc_anot[pbmc_anot.obs['orig.ident']=='MS19270'].obs['labels'].value_counts()
    anot_df
    cluster = anot_df.index.tolist()
    anot_df = anot_df.values/anot_df.values.sum()*100
    anot_df.tolist()
    anot_df = pd.DataFrame({'pct_of_sample':anot_df, 'cluster':cluster})
    anot_df['celltype'] = 'celltype'
    anot_df
    sns.set(rc={'figure.figsize':(8,7)})
    sns.set_style('whitegrid')
    anot_df
    ax = sns.histplot(anot_df, x='celltype', hue='cluster', weights='pct_of_sample', multiple='stack')
    sns.move_legend(ax, loc='center left', bbox_to_anchor=(1.25, 0.5), borderaxespad=0, ncol=2, title=None, frameon=False)
    plt.tight_layout()
    plt.show()

def plot_freq2():
    anot = adata_sub2.obs.copy()
    anot['count'] = 1
    # make variable df for merging later
    variables = pd.DataFrame({'orig.ident':anot['orig.ident']}).drop_duplicates()
    # make df with percentages
    anot_df = anot.groupby(['orig.ident', 'predictions'],sort=False).agg({'count':'sum'}).groupby(level=0).apply(lambda x: 100*x/x.sum()).groupby(level = ['orig.ident', 'predictions']).max().reset_index()
    #merge data
    anot_df = anot_df.merge(variables, left_on = 'orig.ident', right_on = 'orig.ident', how='left').drop_duplicates()
    # new col names
    anot_df = anot_df.set_axis(['orig.ident', 'predictions', 'Pct of sample'], axis=1)
    # plot by condition and celltype
    sns.set_style()
    sns.set(rc={'figure.figsize':(7,6)})
    sns.set_style(rc={'grid.color':'0.9', 'axes.facecolor':'white', 'axes.edgecolor':'grey'})
    plt.ylim(0,100)
    ax = sns.boxplot(x = 'predictions', y='Pct of sample', data =
    anot_df, fliersize = 0, showmeans=True, meanprops={'marker':'+', 'markersize':'10', 'markeredgewidth':'2.3', 'markeredgecolor':'#777777'}, palette = 'turbo')
    ax = sns.stripplot(x= 'predictions', y='Pct of sample', data =
    anot_df, alpha=.5, palette = 'turbo')
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode='anchor')
    ax.tick_params(bottom=True, left=True)
    plt.ylabel('% of sample')
    plt.show()

def get_stats_adata():
    stats_adata_get = simpledialog.askstring(title='reference group',  prompt='What is your reference group to compare conditions, e.g. Control (No quotation marks)')
    stats_adata_var.set(value=stats_adata_get)
    global stats_adata
    stats_adata = stats_adata_var

def make_stats_window(x):
    table_window = Toplevel()
    table_window.title('Results table')
    table_window_res = pt.Table(table_window, dataframe=x, showtoolbar=True, showstatusbar=True)
    table_window_res.showIndex()
    table_window_res.show()

def make_stats():
    print(colored('Calculating stats...', 'yellow'))
    stats_adata_group = stats_adata_var.get()
    lm_res = {}
    celltypes = quant_df.predictions.value_counts().index
    for i in range(0,len(celltypes)):
      df = quant_df[quant_df.predictions==celltypes[i]]
      form = "percentage ~ C(condition, Treatment(reference = '" + stats_adata_group + "'))"
      model = ols(form, data = df)
      fii = model.fit()
      lm_res[i] = fii.summary2().tables[1]
      lm_res[i]['celltype'] = celltypes[i]
      conds = df.condition.value_counts().index.categories.to_list()
      conds
      conds.remove(stats_adata_group)
      conds
      lm_res[i] = lm_res[i].iloc[1:,:]
      comparison = pd.DataFrame(list([stats_adata_group])*len(conds)) + pd.DataFrame(list([' vs '])*len(conds)) + pd.DataFrame(conds)
      comparison[0].values
      lm_res[i]['comparison'] = comparison[0].values
    lm_res2 = pd.concat(lm_res)
    make_stats_window(x=lm_res2)

def plot_freq3():
    # dataframe for plots --> extract values from obs
    anot_df = adata_sub2.obs.copy()
    anot_df['count'] = 1
    # make variable df for merging later
    variables = pd.DataFrame({'orig.ident':anot_df['orig.ident'], 'condition':anot_df['condition']}).drop_duplicates()
    # make df with percentages
    anot_df = anot_df.groupby(['orig.ident', 'predictions', 'condition'],sort=False).agg({'count':'sum'}).groupby(level=0).apply(lambda x: 100*x/x.sum()).groupby(level = ['orig.ident', 'predictions']).max().reset_index()
    #merge data
    anot_df = anot_df.merge(variables, left_on = 'orig.ident', right_on = 'orig.ident', how='left').drop_duplicates()
    # new col names
    anot_df['count'].sum()
    anot_df = anot_df.set_axis(['orig.ident', 'predictions', 'percentage', 'condition'], axis=1)
    anot_df_recall = anot_df.copy() # ! This is to get back to "all subsets" level when aggregation was done for lineages but all subsets need to be plotted afterwards
    # define dataset to perform statistical analysis in stats_function
    global quant_df
    quant_df = anot_df

    # set style parameters
    sns.set_style()
    sns.set(rc={'figure.figsize':(30,7)})
    sns.set_style(rc={'grid.color':'0.9', 'axes.facecolor':'white', 'axes.edgecolor':'grey'})

    # lineages
    if len(anot_df.predictions.value_counts()) > 15:
        #redefine anot_df for lineage
        sub_to_lin = dict({'Th2/Th22':'CD4 T cells', 'Tfh':'CD4 T cells', 'Th1':'CD4 T cells', 'MG CX3CR1+':'Myeloid cells', 'CD8 EM HLA-DRA+':'CD8 T cells', 'CCR5high Th17.1':'CD4 T cells', 'CD4 TEMRA':'CD4 T cells', 'CD8 TRM ITGA1-':'CD8 T cells', 'CD8 CM':'CD8 T cells', 'BAM EMP3+':'Myeloid cells', 'Th17':'CD4 T cells', 'Tregs':'CD4 T cells', 'BAM MRC1+':'Myeloid cells', 'CD8 EM CD160+':'CD8 T cells', 'CD8 TRM ITGA1+':'CD8 T cells', 'mDCs CD1c+': 'Myeloid cells', 'MG CCL2+':'Myeloid cells', 'MG TREM2hi':'Myeloid cells', 'Monocytes':'Myeloid cells', 'CD8 CTL': 'CD8 T cells', 'pDCs':'pDCs', 'NK bright':'NK cells', 'gdT Vd2+':'CD8 T cells',' B activated':'B cells', 'NK dim':'NK cells', 'B atypical':'B cells', 'gdT Vd2-':'CD8 T cells', 'TR-NK':'NK cells', 'MAIT':'CD8 T cells','ILC':'NK cells','Plasmacells':'B cells', 'mDCs CLEC9A+':'Myeloid cells', 'B IL4R+':'B cells', 'mDCs AXL+SIGLEC6+':'Myeloid cells', 'CD8 prolif':'CD8 T cells'})
        adata_sub2.obs['lineage'] = adata_sub2.obs['predictions'].replace(sub_to_lin)
        anot_df = adata_sub2.obs.copy()
        anot_df['count'] = 1
        variables = pd.DataFrame({'orig.ident':anot_df['orig.ident'], 'condition':anot_df['condition']}).drop_duplicates()
        anot_df = anot_df.groupby(['orig.ident', 'lineage', 'condition'],sort=False).agg({'count':'sum'}).groupby(level=0).apply(lambda x: 100*x/x.sum()).groupby(level = ['orig.ident', 'lineage']).max().reset_index()
        anot_df = anot_df.merge(variables, left_on = 'orig.ident', right_on = 'orig.ident', how='left').drop_duplicates()
        anot_df['count'].sum()
        anot_df = anot_df.set_axis(['orig.ident', 'lineage', 'percentage', 'condition'], axis=1)
        # make plot
        ax = sns.boxplot(x = 'lineage', y='percentage', data = anot_df, hue='condition', palette = len(anot_df.condition.value_counts())*["white"], dodge=True, fliersize = 0, showmeans=True, meanprops={'marker':'+', 'markersize':'6', 'markeredgewidth':'1.5', 'markeredgecolor':'#777777'})
        ax = sns.stripplot(x= 'lineage', y='percentage', data = anot_df, alpha=.5, hue = 'condition', dodge=True)
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode='anchor')
        ax.tick_params(bottom=True, left=True)
        plt.ylabel('% of sample')
        plt.show()
        path = filedialog.askdirectory(title = 'Where to store quantification results?') + "/lineage_quantification.csv"
        anot_df.to_csv(path)

    # all subsets
    anot_df = anot_df_recall
    if len(anot_df.predictions.value_counts()) > 15:
        # plot predictions
        ax = sns.boxplot(x = 'predictions', y='percentage', data = anot_df, hue='condition', dodge=True, fliersize = 0, showmeans=True, meanprops={'marker':'+', 'markersize':'6', 'markeredgewidth':'1.5', 'markeredgecolor':'#777777'}, palette = len(anot_df.condition.value_counts())*["white"], order=[' B activated', 'B IL4R+', 'B atypical', 'Plasmacells', 'pDCs', 'Monocytes', 'BAM MRC1+', 'BAM EMP3+', 'MG CX3CR1+', 'MG CCL2+', 'MG TREM2hi', 'mDCs CD1c+', 'mDCs AXL+SIGLEC6+', 'mDCs CLEC9A+', 'NK bright', 'NK dim', 'TR-NK', 'ILC', 'MAIT', 'gdT Vd2+', 'gdT Vd2-', 'CD8 CM', 'CD8 EM HLA-DRA+', 'CD8 EM CD160+', 'CD8 TRM ITGA1+', 'CD8 TRM ITGA1-', 'CD8 CTL', 'Tfh', 'Th17', 'Th2/Th22', 'Th1', 'Tregs', 'CCR5high Th17.1', 'CD4 TEMRA'])
        ax = sns.stripplot(x= 'predictions', y='percentage', data = anot_df, alpha=.5, palette = 'turbo', hue = 'condition', order=[' B activated', 'B IL4R+', 'B atypical', 'Plasmacells', 'pDCs', 'Monocytes', 'BAM MRC1+', 'BAM EMP3+', 'MG CX3CR1+', 'MG CCL2+', 'MG TREM2hi', 'mDCs CD1c+', 'mDCs AXL+SIGLEC6+', 'mDCs CLEC9A+', 'NK bright', 'NK dim', 'TR-NK', 'ILC', 'MAIT', 'gdT Vd2+', 'gdT Vd2-', 'CD8 CM', 'CD8 EM HLA-DRA+', 'CD8 EM CD160+', 'CD8 TRM ITGA1+', 'CD8 TRM ITGA1-', 'CD8 CTL', 'Tfh', 'Th17', 'Th2/Th22', 'Th1', 'Tregs', 'CCR5high Th17.1', 'CD4 TEMRA'], dodge=True)

    # for only subsets:
    if 'CD8 CM' in anot_df.predictions.values and 'Th2/Th22' not in anot_df.predictions.values:
        ax = sns.boxplot(x = 'predictions', y='percentage', data = anot_df, hue='condition', dodge=True, fliersize = 0, showmeans=True, meanprops={'marker':'+', 'markersize':'6', 'markeredgewidth':'1.5', 'markeredgecolor':'#777777'}, palette = 'turbo', order=['MAIT', 'gdT Vd2+', 'gdT Vd2-', 'CD8 CM', 'CD8 EM HLA-DRA+', 'CD8 EM CD160+', 'CD8 TRM ITGA1+', 'CD8 TRM ITGA1-', 'CD8 CTL'])
        ax = sns.stripplot(x= 'predictions', y='percentage', data = anot_df, alpha=.5, palette = len(anot_df.condition.value_counts())*["white"], hue = 'condition', dodge=True, order= ['MAIT', 'gdT Vd2+', 'gdT Vd2-', 'CD8 CM', 'CD8 EM HLA-DRA+', 'CD8 EM CD160+', 'CD8 TRM ITGA1+', 'CD8 TRM ITGA1-', 'CD8 CTL'])
    elif ' B activated' in anot_df.predictions.values and 'Th2/Th22' not in anot_df.predictions.values:
        ax = sns.boxplot(x = 'predictions', y='percentage', data = anot_df, hue='condition', dodge=True, fliersize = 0, showmeans=True, meanprops={'marker':'+', 'markersize':'6', 'markeredgewidth':'1.5', 'markeredgecolor':'#777777'}, palette = 'turbo', order=[' B activated', 'B IL4R+', 'B atypical', 'Plasmacells'])
        ax = sns.stripplot(x= 'predictions', y='percentage', data = anot_df, alpha=.5, palette = len(anot_df.condition.value_counts())*["white"], hue = 'condition', dodge=True, order = [' B activated', 'B IL4R+', 'B atypical', 'Plasmacells'])
    elif 'CD4 TEMRA' in anot_df.predictions.values and 'BAM EMP3+' not in anot_df.predictions.values:
        ax = sns.boxplot(x = 'predictions', y='percentage', data = anot_df, hue='condition', dodge=True, fliersize = 0, showmeans=True, meanprops={'marker':'+', 'markersize':'6', 'markeredgewidth':'1.5', 'markeredgecolor':'#777777'}, palette = 'turbo', order=['Tfh', 'Th17', 'Th2/Th22', 'Th1', 'Tregs', 'CCR5high Th17.1', 'CD4 TEMRA'])
        ax = sns.stripplot(x= 'predictions', y='percentage', data = anot_df, alpha=.5, palette = len(anot_df.condition.value_counts())*["white"], hue = 'condition', dodge=True, order=['Tfh', 'Th17', 'Th2/Th22', 'Th1', 'Tregs', 'CCR5high Th17.1', 'CD4 TEMRA'])
    elif 'BAM MRC1+' in anot_df.predictions.values and 'Th2/Th22' not in anot_df.predictions.values:
        ax = sns.boxplot(x = 'predictions', y='percentage', data = anot_df, hue='condition', dodge=True, fliersize = 0, showmeans=True, meanprops={'marker':'+', 'markersize':'6', 'markeredgewidth':'1.5', 'markeredgecolor':'#777777'}, palette = 'turbo', order=['Monocytes', 'BAM MRC1+', 'BAM EMP3+', 'MG CX3CR1+', 'MG CCL2+', 'MG TREM2hi', 'mDCs CD1c+', 'mDCs AXL+SIGLEC6+', 'mDCs CLEC9A+'])
        ax = sns.stripplot(x= 'predictions', y='percentage', data = anot_df, alpha=.5, palette = len(anot_df.condition.value_counts())*["white"], hue = 'condition', dodge=True, order=['Monocytes', 'BAM MRC1+', 'BAM EMP3+', 'MG CX3CR1+', 'MG CCL2+', 'MG TREM2hi', 'mDCs CD1c+', 'mDCs AXL+SIGLEC6+', 'mDCs CLEC9A+'])
    elif 'NK bright' in anot_df.predictions.values and 'Th2/Th22' not in anot_df.predictions.values:
        ax = sns.boxplot(x = 'predictions', y='percentage', data = anot_df, hue='condition', dodge=True, fliersize = 0, showmeans=True, meanprops={'marker':'+', 'markersize':'6', 'markeredgewidth':'1.5', 'markeredgecolor':'#777777'}, palette = 'turbo', order=['NK bright', 'NK dim', 'TR-NK', 'ILC'])
        ax = sns.stripplot(x= 'predictions', y='percentage', data = anot_df, alpha=.5, palette = len(anot_df.condition.value_counts())*["white"], hue = 'condition', dodge=True, order=['NK bright', 'NK dim', 'TR-NK', 'ILC'])

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode='anchor')
    ax.tick_params(bottom=True, left=True)
    plt.ylabel('% of sample')
    plt.show()
    path = filedialog.askdirectory(title = 'Where to store quantification results?') + "/subset_quantification.csv"
    anot_df.to_csv(path)

    # ask whether to change color or not
    print(colored('Close the plot to proceed. You can now choose to change default colors and order of conditions. To view a subset of conditions (e.g. only Control), only type the name of this condition and a color of your choice.', 'yellow'))
    answer = messagebox.askyesno(title='Change plot', message='Do you want to change colors and condition order?')

    if answer == True:
        # get inputted category order
        get_condition_order()
        # get inputter hue order
        get_hues()
        # read new condition order
        new_order = condition_order.get()
        new_order = eval(new_order)
        new_order = list(new_order)
        # read new colors
        hues = new_hues.get()
        hues = eval(hues)
        colors = sns.set_palette(sns.color_palette(hues))
        # if no subset was chosen, also show plot by lineage
        anot_df = anot_df_recall
        if len(anot_df.predictions.value_counts()) > 15:
            #redefine anot_df for lineage
            sub_to_lin = dict({'Th2/Th22':'CD4 T cells', 'Tfh':'CD4 T cells', 'Th1':'CD4 T cells', 'MG CX3CR1+':'Myeloid cells', 'CD8 EM HLA-DRA+':'CD8 T cells', 'CCR5high Th17.1':'CD4 T cells', 'CD4 TEMRA':'CD4 T cells', 'CD8 TRM ITGA1-':'CD8 T cells', 'CD8 CM':'CD8 T cells', 'BAM EMP3+':'Myeloid cells', 'Th17':'CD4 T cells', 'Tregs':'CD4 T cells', 'BAM MRC1+':'Myeloid cells', 'CD8 EM CD160+':'CD8 T cells', 'CD8 TRM ITGA1+':'CD8 T cells', 'mDCs CD1c+': 'Myeloid cells', 'MG CCL2+':'Myeloid cells', 'MG TREM2hi':'Myeloid cells', 'Monocytes':'Myeloid cells', 'CD8 CTL': 'CD8 T cells', 'pDCs':'pDCs', 'NK bright':'NK cells', 'gdT Vd2+':'CD8 T cells',' B activated':'B cells', 'NK dim':'NK cells', 'B atypical':'B cells', 'gdT Vd2-':'CD8 T cells', 'TR-NK':'NK cells', 'MAIT':'CD8 T cells','ILC':'NK cells','Plasmacells':'B cells', 'mDCs CLEC9A+':'Myeloid cells', 'B IL4R+':'B cells', 'mDCs AXL+SIGLEC6+':'Myeloid cells', 'CD8 prolif':'CD8 T cells'})
            adata_sub2.obs['lineage'] = adata_sub2.obs['predictions'].replace(sub_to_lin)
            anot_df = adata_sub2.obs.copy()
            anot_df['count'] = 1
            variables = pd.DataFrame({'orig.ident':anot_df['orig.ident'], 'condition':anot_df['condition']}).drop_duplicates()
            anot_df = anot_df.groupby(['orig.ident', 'lineage', 'condition'],sort=False).agg({'count':'sum'}).groupby(level=0).apply(lambda x: 100*x/x.sum()).groupby(level = ['orig.ident', 'lineage']).max().reset_index()
            anot_df = anot_df.merge(variables, left_on = 'orig.ident', right_on = 'orig.ident', how='left').drop_duplicates()
            anot_df['count'].sum()
            anot_df = anot_df.set_axis(['orig.ident', 'lineage', 'percentage', 'condition'], axis=1)
            # make plot
            ax = sns.boxplot(x = 'lineage', y='percentage', data = anot_df, hue='condition',  dodge=True, fliersize = 0, showmeans=True, hue_order=new_order, palette= len(anot_df.condition.value_counts())*["white"], meanprops={'marker':'+', 'markersize':'6', 'markeredgewidth':'1.5', 'markeredgecolor':'#777777'})
            ax = sns.stripplot(x= 'lineage', y='percentage', data = anot_df, alpha=.5, hue = 'condition', dodge=True, hue_order=new_order, palette=colors)
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode='anchor')
            ax.tick_params(bottom=True, left=True)
            plt.ylabel('% of sample')
            plt.show()
        # plot only necessary subsets
        anot_df = anot_df_recall
        if len(anot_df.predictions.value_counts()) > 15:
            anot_df = anot_df_recall
            ax = sns.boxplot(x = 'predictions', y='percentage', data = anot_df, hue='condition', hue_order = new_order,  dodge=True, fliersize = 0, showmeans=True, meanprops={'marker':'+', 'markersize':'6', 'markeredgewidth':'1.5', 'markeredgecolor':'#777777'}, palette = len(anot_df.condition.value_counts())*["white"], order=[' B activated', 'B IL4R+', 'B atypical', 'Plasmacells', 'pDCs', 'Monocytes', 'BAM MRC1+', 'BAM EMP3+', 'MG CX3CR1+', 'MG CCL2+', 'MG TREM2hi', 'mDCs CD1c+', 'mDCs AXL+SIGLEC6+', 'mDCs CLEC9A+', 'NK bright', 'NK dim', 'TR-NK', 'ILC', 'MAIT', 'gdT Vd2+', 'gdT Vd2-', 'CD8 CM', 'CD8 EM HLA-DRA+', 'CD8 EM CD160+', 'CD8 TRM ITGA1+', 'CD8 TRM ITGA1-', 'CD8 CTL', 'Tfh', 'Th17', 'Th2/Th22', 'Th1', 'Tregs', 'CCR5high Th17.1', 'CD4 TEMRA'])
            ax = sns.stripplot(x= 'predictions', y='percentage', data = anot_df, alpha=.5, palette = colors, hue = 'condition', order=[' B activated', 'B IL4R+', 'B atypical', 'Plasmacells', 'pDCs', 'Monocytes', 'BAM MRC1+', 'BAM EMP3+', 'MG CX3CR1+', 'MG CCL2+', 'MG TREM2hi', 'mDCs CD1c+', 'mDCs AXL+SIGLEC6+', 'mDCs CLEC9A+', 'NK bright', 'NK dim', 'TR-NK', 'ILC', 'MAIT', 'gdT Vd2+', 'gdT Vd2-', 'CD8 CM', 'CD8 EM HLA-DRA+', 'CD8 EM CD160+', 'CD8 TRM ITGA1+', 'CD8 TRM ITGA1-', 'CD8 CTL', 'Tfh', 'Th17', 'Th2/Th22', 'Th1', 'Tregs', 'CCR5high Th17.1', 'CD4 TEMRA'], dodge=True, hue_order = new_order)

        # for only subsets:
        if 'CD8 CM' in anot_df.predictions.values and 'Th2/Th22' not in anot_df.predictions.values:
            ax = sns.boxplot(x = 'predictions', y='percentage', data = anot_df, hue='condition', hue_order = new_order, dodge=True, fliersize = 0, showmeans=True, meanprops={'marker':'+', 'markersize':'6', 'markeredgewidth':'1.5', 'markeredgecolor':'#777777'}, palette = len(anot_df.condition.value_counts())*["white"], order=['MAIT', 'gdT Vd2+', 'gdT Vd2-', 'CD8 CM', 'CD8 EM HLA-DRA+', 'CD8 EM CD160+', 'CD8 TRM ITGA1+', 'CD8 TRM ITGA1-', 'CD8 CTL'])
            ax = sns.stripplot(x= 'predictions', y='percentage', data = anot_df, alpha=.5, palette = colors, hue_order = new_order, hue = 'condition', dodge=True, order=['MAIT', 'gdT Vd2+', 'gdT Vd2-', 'CD8 CM', 'CD8 EM HLA-DRA+', 'CD8 EM CD160+', 'CD8 TRM ITGA1+', 'CD8 TRM ITGA1-', 'CD8 CTL'])
        elif ' B activated' in anot_df.predictions.values and 'Th2/Th22' not in anot_df.predictions.values:
            ax = sns.boxplot(x = 'predictions', y='percentage', data = anot_df, hue='condition', hue_order = new_order, dodge=True, fliersize = 0, showmeans=True, meanprops={'marker':'+', 'markersize':'6', 'markeredgewidth':'1.5', 'markeredgecolor':'#777777'}, palette = len(anot_df.condition.value_counts())*["white"], order=[' B activated', 'B IL4R+', 'B atypical', 'Plasmacells'])
            ax = sns.stripplot(x= 'predictions', y='percentage', data = anot_df, alpha=.5, palette = colors, hue_order = new_order, hue = 'condition', dodge=True, order=[' B activated', 'B IL4R+', 'B atypical', 'Plasmacells'])
        elif 'CD4 TEMRA' in anot_df.predictions.values and 'BAM EMP3+' not in anot_df.predictions.values:
            ax = sns.boxplot(x = 'predictions', y='percentage', data = anot_df, hue='condition', hue_order = new_order, dodge=True, fliersize = 0, showmeans=True, meanprops={'marker':'+', 'markersize':'6', 'markeredgewidth':'1.5', 'markeredgecolor':'#777777'}, palette = len(anot_df.condition.value_counts())*["white"], order=['Tfh', 'Th17', 'Th2/Th22', 'Th1', 'Tregs', 'CCR5high Th17.1', 'CD4 TEMRA'])
            ax = sns.stripplot(x= 'predictions', y='percentage', data = anot_df, alpha=.5, palette = colors, hue_order = new_order, hue = 'condition', dodge=True, order=['Tfh', 'Th17', 'Th2/Th22', 'Th1', 'Tregs', 'CCR5high Th17.1', 'CD4 TEMRA'])
        elif 'BAM MRC1+' in anot_df.predictions.values and 'Th2/Th22' not in anot_df.predictions.values:
            ax = sns.boxplot(x = 'predictions', y='percentage', data = anot_df, hue='condition', hue_order = new_order, dodge=True, fliersize = 0, showmeans=True, meanprops={'marker':'+', 'markersize':'6', 'markeredgewidth':'1.5', 'markeredgecolor':'#777777'}, palette = len(anot_df.condition.value_counts())*["white"], order=['Monocytes', 'BAM MRC1+', 'BAM EMP3+', 'MG CX3CR1+', 'MG CCL2+', 'MG TREM2hi', 'mDCs CD1c+', 'mDCs AXL+SIGLEC6+', 'mDCs CLEC9A+'])
            ax = sns.stripplot(x= 'predictions', y='percentage', data = anot_df, alpha=.5, palette = colors, hue_order = new_order, hue = 'condition', dodge=True, order=['Monocytes', 'BAM MRC1+', 'BAM EMP3+', 'MG CX3CR1+', 'MG CCL2+', 'MG TREM2hi', 'mDCs CD1c+', 'mDCs AXL+SIGLEC6+', 'mDCs CLEC9A+'])
        elif 'NK bright' in anot_df.predictions.values and 'Th2/Th22' not in anot_df.predictions.values:
            ax = sns.boxplot(x = 'predictions', y='percentage', data = anot_df, hue='condition', hue_order = new_order, dodge=True, fliersize = 0, showmeans=True, meanprops={'marker':'+', 'markersize':'6', 'markeredgewidth':'1.5', 'markeredgecolor':'#777777'}, palette = len(anot_df.condition.value_counts())*["white"], order=['NK bright', 'NK dim', 'TR-NK', 'ILC'])
            ax = sns.stripplot(x= 'predictions', y='percentage', data = anot_df, alpha=.5, palette = colors, hue_order = new_order, hue = 'condition', dodge=True, order=['NK bright', 'NK dim', 'TR-NK', 'ILC'])

        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode='anchor')
        ax.tick_params(bottom=True, left=True)
        plt.ylabel('% of sample')
        plt.show()

        # if no subset was chosen, also show plot by lineage
        anot_df = anot_df_recall
        if len(anot_df.predictions.value_counts()) > 15:
            #redefine anot_df for lineage
            sub_to_lin = dict({'Th2/Th22':'CD4 T cells', 'Tfh':'CD4 T cells', 'Th1':'CD4 T cells', 'MG CX3CR1+':'Myeloid cells', 'CD8 EM HLA-DRA+':'CD8 T cells', 'CCR5high Th17.1':'CD4 T cells', 'CD4 TEMRA':'CD4 T cells', 'CD8 TRM ITGA1-':'CD8 T cells', 'CD8 CM':'CD8 T cells', 'BAM EMP3+':'Myeloid cells', 'Th17':'CD4 T cells', 'Tregs':'CD4 T cells', 'BAM MRC1+':'Myeloid cells', 'CD8 EM CD160+':'CD8 T cells', 'CD8 TRM ITGA1+':'CD8 T cells', 'mDCs CD1c+': 'Myeloid cells', 'MG CCL2+':'Myeloid cells', 'MG TREM2hi':'Myeloid cells', 'Monocytes':'Myeloid cells', 'CD8 CTL': 'CD8 T cells', 'pDCs':'pDCs', 'NK bright':'NK cells', 'gdT Vd2+':'CD8 T cells',' B activated':'B cells', 'NK dim':'NK cells', 'B atypical':'B cells', 'gdT Vd2-':'CD8 T cells', 'TR-NK':'NK cells', 'MAIT':'CD8 T cells','ILC':'NK cells','Plasmacells':'B cells', 'mDCs CLEC9A+':'Myeloid cells', 'B IL4R+':'B cells', 'mDCs AXL+SIGLEC6+':'Myeloid cells', 'CD8 prolif':'CD8 T cells'})
            adata_sub2.obs['lineage'] = adata_sub2.obs['predictions'].replace(sub_to_lin)
            anot_df = adata_sub2.obs.copy()
            anot_df['count'] = 1
            variables = pd.DataFrame({'orig.ident':anot_df['orig.ident'], 'condition':anot_df['condition']}).drop_duplicates()
            anot_df = anot_df.groupby(['orig.ident', 'lineage', 'condition'],sort=False).agg({'count':'sum'}).groupby(level=0).apply(lambda x: 100*x/x.sum()).groupby(level = ['orig.ident', 'lineage']).max().reset_index()
            anot_df = anot_df.merge(variables, left_on = 'orig.ident', right_on = 'orig.ident', how='left').drop_duplicates()
            anot_df['count'].sum()
            anot_df = anot_df.set_axis(['orig.ident', 'lineage', 'percentage', 'condition'], axis=1)
            # make plot
            ax = sns.boxplot(x = 'lineage', y='percentage', data = anot_df, hue='condition', hue_order=hues, dodge=True, fliersize = 0, showmeans=True, meanprops={'marker':'+', 'markersize':'6', 'markeredgewidth':'1.5', 'markeredgecolor':'#777777'}, palette = len(anot_df.condition.value_counts())*["white"])
            ax = sns.stripplot(x= 'lineage', y='percentage', data = anot_df, alpha=.5, palette = colors, hue = 'condition', dodge=True, hue_order=hues)
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode='anchor')
            ax.tick_params(bottom=True, left=True)
            plt.ylabel('% of sample')
            plt.show()

        get_stats_adata()
        make_stats()
    else:
        get_stats_adata()
        make_stats()

def write_data():
    data_type = data_to_write.get()
    get_path()
    if data_type == '.h5ad':
        print(colored('Writing to .h5ad...', 'red'))
        #adata_to_sub.var.index.name = "gene symbols"
        from scipy.sparse import csr_matrix
        adata_to_sub.X = csr_matrix(adata_to_sub.X)
        path = path_var.get()
        path = path + "/labelled_query.h5ad"
        print('Writing data to:' + path)
        adata_to_sub.write_h5ad(path)
        print(colored('...done!', 'blue'))
    elif data_type == 'Annotations':
        print(colored('Writing annotations to .tsv...', 'red'))
        csv_out = pd.DataFrame({'mapped_labels':adata_to_sub.obs.predictions, 'orig.ident':adata_to_sub.obs['orig.ident'], 'study':adata_to_sub.obs['study'], 'condition':adata_to_sub.obs['condition']})
        path = path_var.get()
        path = path + "/mapped_labels.tsv"
        print('Writing data to:' + path)
        csv_out.to_csv(path, sep='\t')
        print(colored('...done!', 'cyan'))
    elif data_type == 'UMAP coordinates':
        print(colored('Writing umap coords to .tsv...', 'red'))
        umap_out = pd.DataFrame(adata_to_sub.obsm['X_umap'])
        umap_out = umap_out.rename({0: 'UMAP_1', 1: 'UMAP_2'}, axis=1)
        path = path_var.get()
        path = path + "/umap_coords.tsv"
        umap_out.to_csv(path, sep='\t')
        print(colored('...done!', 'cyan'))

# main scvi workflow

def load_user_adata():

    print(colored('Loading 10X data...', 'red'))
    # load mtx file
    matrix_path = path_var.get() + '/matrix.mtx'
    adata = sc.read(matrix_path).T
    print(colored('...done!', 'cyan'))

    # load gene ids and symbols for annotation
    #check whether features or genes.tsv is provided
    files = os.listdir(path_var.get())
    if 'features' in str(files):
      genes_path = path_var.get() + '/features.tsv'
    else:
      genes_path = path_var.get() + '/genes.tsv'

    # read genes
    genes = pd.read_csv(genes_path, header=None, sep='\t')
    adata.var['gene_ids'] = genes[0].values
    adata.var['gene_symbols'] = genes[1].values
    adata.var_names = adata.var['gene_symbols']
    adata.var_names_make_unique(join="-")

    # load barcodes
    barcodes_path = path_var.get() + '/barcodes.tsv'
    cells = pd.read_csv(barcodes_path, header=None, sep='\t')
    adata.obs['barcode'] = cells[0].values
    adata.obs_names = cells[0]
    adata.obs_names_make_unique(join="-")

    # assign tissue to csf
    adata.obs['tissue'] = 'csf'
    # remove ribosomal protein genes
    adata.var['rp'] = adata.var_names.str.startswith('RP')
    # exclude all ribosomal protein genes
    adata = adata[:,adata.var['rp']==False]
    # make raw count layer
    adata.layers["counts"] = adata.X.copy()

    # add metadata
    if 'idents' in globals():
        adata.obs['orig.ident'] = idents['x'].values
    else:
    	adata.obs['orig.ident'] = "sample"

    if 'condition' in globals():
    	adata.obs['condition'] = condition['x'].values
    else:
    	adata.obs['condition'] = "condition"

    if 'study_ident' in globals():
    	adata.obs['study'] = study_ident['x'].values
    else:
    	adata.obs['study'] = "study"

    # find doublets
    print(colored('Finding doublets...', 'red'))
    scrub_tbl = scr.Scrublet(adata.layers["counts"])
    doublet_scores, predicted_doublets = scrub_tbl.scrub_doublets()
    #sns.histplot(doublet_scores)
    adata.obs['doublets'] = predicted_doublets
    adata.obs['doublet_score'] = doublet_scores
    adata = adata[adata.obs.doublet_score < 0.2]
    # batch_list = []
    # # make vector of sample names
    # samples = adata.obs["orig.ident"].unique()
    # for x in samples:
    #     print('starting', x)
    #     batch = adata[adata.obs["orig.ident"].isin([x])]
    #     scrub = scr.Scrublet(batch.layers["counts"])
    #     doublet_scores, predicted_doublets = scrub.scrub_doublets()
    #     print('doublets found')
    #     batch.obs['doublets'] = predicted_doublets
    #     batch.obs['doublet_score'] = doublet_scores
    #     print('added to object')
    #     batch_list.append(batch)
    #     print(colored('Batch', 'cyan'), colored(x, 'white'), colored(' done!', 'cyan'))
    # # merge datasets
    # adata = batch_list[0]
    # for x in range(1, (len(batch_list))):
    #     adata = adata.concatenate(batch_list[x])
    # adata = adata[adata.obs.doublet_score < 0.2]
    print(colored('Doublets found!', 'cyan'))

    # normalization for HVG calculation
    print(colored('Preprocessing and QC...', 'red'))
    print('DEBUG: normalization 1')
    sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed = True)
    print('DEBUG: normalization 2')
    sc.pp.log1p(adata)
    print('DEBUG: normalization 3')
    adata.raw = adata
    print('DEBUG: normalization 4')

    # QC vars
    print('DEBUG: QC 1')
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    print('DEBUG: QC 2')
    adata.var['rp'] = adata.var_names.str.startswith('RPS', 'RPL')
    print('DEBUG: QC 3')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'rp'], percent_top=None, inplace=True)
    print('DEBUG: QC 4')

    # ask thresholds
    print(colored('Type threshold for QC variables...', 'yellow'))
    ask_threshs()

    # make QC plots
    print('DEBUG: Make QC plots...')
    plt.rcParams['figure.figsize'] = 8,8
    fig, axes = plt.subplots(2,2)
    sns.violinplot(data=adata.obs, x="condition", y="total_counts", ax=axes[0,0])
    sns.violinplot(data=adata.obs, x="condition", y="log1p_total_counts", ax=axes[0,1])
    sns.violinplot(data=adata.obs, x="condition", y="n_genes_by_counts", ax=axes[1,0])
    sns.violinplot(data=adata.obs, x="condition", y="pct_counts_mt", ax=axes[1,1])
    plt.tight_layout()
    print('DEBUG: Make QC plots SHOW')
    plt.show()
    print('DEBUG: Make QC plots finished')

    #sc.pl.violin(adata, keys=['log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_mt'])
    adata = adata[(adata.obs.total_counts > count_thresh.get()) & (adata.obs.n_genes_by_counts > feat_thresh.get()) & (adata.obs.pct_counts_mt < mt_thresh.get())]
    #sc.pl.violin(adata, keys=['log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_mt'])
    print(colored('...done!', 'cyan'))

    # load reference
    print(colored('Loading reference data...', 'red'))
    if getattr(sys, 'frozen', False):
        ref = sc.read_h5ad(os.path.join(sys._MEIPASS, "files/ref2.h5ad"))
    else:
        ref = sc.read_h5ad("files/ref2.h5ad")
    print(colored('Reference loaded!', 'cyan'))

    # assign batch
    ref.obs['batch'] = '0'
    adata.obs['batch'] = '1'

    # rename CD8 TRM labels
    ref.obs['labels'] = np.where(ref.obs['labels']=='CD8 TRM1', 'CD8 TRM ITGA1+', ref.obs['labels'])
    ref.obs['labels'] = np.where(ref.obs['labels']=='CD8 TRM2', 'CD8 TRM ITGA1-', ref.obs['labels'])

    # "ref and query data loaded"
    #message_info2()
    #print(ref.obs['orig.ident'].value_counts())

    # preprocess referece
    sc.pp.normalize_total(ref, target_sum=1e4, exclude_highly_expressed = True)
    sc.pp.log1p(ref)
    ref.raw = ref

    # combine, find hvgs, and split
    print(colored('Concatenating datasets and finding variable genes...', 'red'))
    combined = ref.concatenate(adata)
    sc.pp.highly_variable_genes(combined, n_top_genes=5000, subset=True, flavor = 'seurat', batch_key='orig.ident')
    print('combined dataset:')
    print(combined)
    ref_red = combined[combined.obs.batch=='0'].copy()
    adata_red = combined[combined.obs.batch=='1'].copy()
    print(colored('...done!', 'cyan'))

    # compute qc variables
    ref_red.var['mt'] = ref_red.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(ref_red, qc_vars=['mt'], percent_top=None, inplace=True)
    adata_red.var['mt'] = adata_red.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata_red, qc_vars=['mt'], percent_top=None, inplace=True)

    # define number of hidden neurons and latent variables
    neurons = [160]
    latents = [30]
    print(colored('Setting hyperparameters for VAE. Number of neurons per layer:', 'white'), colored(neurons, 'yellow'), colored('Number of latent variables:', 'white'), colored(latents, 'yellow'))

    # set seed
    random.seed()
    scvi.settings.seed = random.randint(0,99999999)

    # make scVI VAE for ref
    print(colored('Preparing data and commencing training on reference dataset', 'red'))
    scvi.data.setup_anndata(ref_red, layer="counts", batch_key="orig.ident", categorical_covariate_keys = ['tissue', 'study'], continuous_covariate_keys = ['pct_counts_mt'], labels_key='labels')
    ref_model = scvi.model.SCVI(ref_red, n_hidden = neurons[0], n_layers = 2, n_latent = latents[0], dropout_rate=.2)

    # "Start training"
    #message_info3()

    ref_model.train(early_stopping=True, max_epochs = 20)
    print(colored('Training complete.', 'cyan'))
    # train_elbo = ref_model.history['elbo_train'][1:]
    # test_elbo = ref_model.history['elbo_validation']
    # ax = train_elbo.plot()
    # test_elbo.plot(ax=ax)
    # plt.rcParams['figure.figsize'] = 5,5
    # plt.show()
    ref_latent = ref_model.get_latent_representation()
    ref_red.obsm['X_scVI'] = ref_latent
    sc.pp.neighbors(ref_red, use_rep='X_scVI')
    sc.tl.umap(ref_red)
    sns.set(rc={'figure.figsize':(11,11)})
    sns.set_style('whitegrid')
    #make scANVI VAE for ref
    print(colored('Supervising trained model using labels...', 'red'))
    vae_scanvi_model = scvi.model.SCANVI.from_scvi_model(ref_model, unlabeled_category='Unknown', adata=ref_red)

    # "Supervise reference"
    #message_info4()

    vae_scanvi_model.train(max_epochs = 50, early_stopping = True)
    print(colored('...done!', 'cyan'))
    ref_red.obsm['X_scANVI'] = vae_scanvi_model.get_latent_representation()
    sc.pp.neighbors(ref_red, use_rep='X_scANVI')
    sc.tl.umap(ref_red)
    sns.set(rc={'figure.figsize':(11,11)})
    sns.set_style('whitegrid')
    #sc.pl.umap(ref_red, color=['labels'], size=40, legend_loc="on data")

    # map labels to query adata_red
    adata_red.obs['labels'] = 'Unknown'
    ref_red.obs['predictions'] = vae_scanvi_model.predict()
    #sc.pl.umap(ref_red, color=['predictions'], size=40, legend_loc="on data")

    global vae_adata_red
    vae_adata_red = scvi.model.SCANVI.load_query_data(adata_red, vae_scanvi_model)

    # "Mapping labels"
    #message_info5()
    print(colored('Starting mapping process... This may take a while. Using GPU via cuda recommended.', 'red'))
    vae_adata_red.train(max_epochs = 400, early_stopping = True, check_val_every_n_epoch = 1, plan_kwargs=dict(weight_decay=0.0), early_stopping_min_delta =0.01, early_stopping_patience = 15)
    adata_red.obsm['X_scANVI'] = vae_adata_red.get_latent_representation()
    adata_red.obs['predictions'] = vae_adata_red.predict()
    adata_red.obs.predictions
    sc.pp.neighbors(adata_red, use_rep = 'X_scANVI')
    sc.tl.umap(adata_red)
    #check in query
    sns.set(rc={'figure.figsize':(15,15)})
    sns.set_style('whitegrid')

    # "Mapping complete!"
    print(colored('-------------------------------------------------------------------------\n', 'white'), colored('Mapping complete! You can now continue with plotting and quantification.', 'green'))
    # declare object as global for use in plotting and writing functions, also make copy to recall after subetting
    global adata_to_sub
    global bdata
    adata_to_sub = adata_red.copy()
    bdata = adata_red.copy()
    bdata.layers['scvi_norm'] = vae_adata_red.get_normalized_expression(library_size=10e4)

def load_user_adata2():
    print(colored('Loading 10X data...', 'red'))
    # load mtx file
    matrix_path = path_var.get() + '/matrix.mtx'
    adata = sc.read(matrix_path).T
    print(colored('...done!', 'cyan'))
    # load gene ids and symbols for annotation
    #check whether features or genes.tsv is provided
    files = os.listdir(path_var.get())
    if 'features' in str(files):
      genes_path = path_var.get() + '/features.tsv'
    else:
      genes_path = path_var.get() + '/genes.tsv'
    # read genes
    genes = pd.read_csv(genes_path, header=None, sep='\t')
    adata.var['gene_ids'] = genes[0].values
    adata.var['gene_symbols'] = genes[1].values
    adata.var_names = adata.var['gene_symbols']
    adata.var_names_make_unique(join="-")
    # load barcodes
    barcodes_path = path_var.get() + '/barcodes.tsv'
    cells = pd.read_csv(barcodes_path, header=None, sep='\t')
    adata.obs['barcode'] = cells[0].values
    adata.obs_names = cells[0]
    adata.obs_names_make_unique(join="-")
    # assign tissue to csf
    adata.obs['tissue'] = 'csf'
    # remove ribosomal protein genes
    adata.var['rp'] = adata.var_names.str.startswith('RP')
    # exclude all ribosomal protein genes
    adata = adata[:,adata.var['rp']==False]
    # make raw count layer
    adata.layers["counts"] = adata.X.copy()
    # add metadata
    if 'idents' in globals():
        adata.obs['orig.ident'] = idents['x'].values
    else:
    	adata.obs['orig.ident'] = "sample"

    if 'condition' in globals():
    	adata.obs['condition'] = condition['x'].values
    else:
    	adata.obs['condition'] = "condition"

    if 'study_ident' in globals():
    	adata.obs['study'] = study_ident['x'].values
    else:
    	adata.obs['study'] = "study"
    # find doublets
    print(colored('Finding doublets...', 'red'))
    scrub_tbl = scr.Scrublet(adata.layers["counts"])
    doublet_scores, predicted_doublets = scrub_tbl.scrub_doublets()
    #sns.histplot(doublet_scores)
    adata.obs['doublets'] = predicted_doublets
    adata.obs['doublet_score'] = doublet_scores
    adata = adata[adata.obs.doublet_score < 0.2]
    # batch_list = []
    # # make vector of sample names
    # samples = adata.obs["orig.ident"].unique()
    # for x in samples:
    #     print('starting', x)
    #     batch = adata[adata.obs["orig.ident"].isin([x])]
    #     scrub = scr.Scrublet(batch.layers["counts"])
    #     doublet_scores, predicted_doublets = scrub.scrub_doublets()
    #     print('doublets found')
    #     batch.obs['doublets'] = predicted_doublets
    #     batch.obs['doublet_score'] = doublet_scores
    #     print('added to object')
    #     batch_list.append(batch)
    #     print(colored('Batch', 'cyan'), colored(x, 'white'), colored(' done!', 'cyan'))
    # # merge datasets
    # adata = batch_list[0]
    # for x in range(1, (len(batch_list))):
    #     adata = adata.concatenate(batch_list[x])
    # adata = adata[adata.obs.doublet_score < 0.2]
    print(colored('Doublets found!', 'cyan'))
    # normalization for HVG calculation
    print(colored('Preprocessing and QC...', 'red'))
    # QC vars
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['rp'] = adata.var_names.str.startswith('RPS', 'RPL')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'rp'], percent_top=None, inplace=True)
    # ask hyperolds
    print(adata.obs.condition.value_counts())
    print(adata.obs.tissue.value_counts())
    print(adata.obs['orig.ident'].value_counts())
    print(adata.obs)
    ask_threshs()
    # make QC plots
    plt.rcParams['figure.figsize'] = 8,8
    fig, axes = plt.subplots(2,2)
    sns.violinplot(data=adata.obs, x="condition", y="total_counts", ax=axes[0,0])
    sns.violinplot(data=adata.obs, x="condition", y="n_genes_by_counts", ax=axes[0,1])
    sns.violinplot(data=adata.obs, x="condition", y="pct_counts_mt", ax=axes[1,0])
    sns.violinplot(data=adata.obs, x="condition", y="pct_counts_mt", ax=axes[1,1])
    plt.show()
    # # get qc hyperolds from user
    #sc.pl.violin(adata, keys=['log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_mt'])
    print(adata)
    adata = adata[(adata.obs.total_counts > count_thresh.get()) & (adata.obs.n_genes_by_counts > feat_thresh.get()) & (adata.obs.pct_counts_mt < mt_thresh.get())]
    print(adata)
    #sc.pl.violin(adata, keys=['log1p_n_genes_by_counts', 'log1p_total_counts', 'pct_counts_mt'])
    print(colored('...done!', 'cyan'))
    # preprocess adataerece
    sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed = True)
    sc.pp.log1p(adata)
    adata.raw = adata
    # combine, find hvgs, and split
    print(colored('Finding variable genes...', 'red'))
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=True, flavor = 'seurat', batch_key='orig.ident')
    # compute qc variables
    # define number of hidden neurons and latent variables
    ask_hypers()
    neurons = num_neurons.get()
    layers = num_layers.get()
    latents = num_latents.get()
    print(colored('Setting hyperparameters for VAE. Number of neurons per layer:', 'white'), colored(neurons, 'yellow'), colored('Number of layers:', 'white'), colored(layers, 'yellow'), colored('Number of latent variables:', 'white'), colored(latents, 'yellow'))
    # set seed
    random.seed()
    scvi.settings.seed = random.randint(0,99999999)
    # make scVI VAE for adata
    print(colored('Preparing data and commencing training on dataset', 'red'))
    scvi.data.setup_anndata(adata, layer="counts", batch_key="orig.ident", categorical_covariate_keys = ['tissue', 'study', 'condition'], continuous_covariate_keys = ['pct_counts_mt'])
    global vae_adata_red
    vae_adata_red = scvi.model.SCVI(adata, n_hidden = neurons, n_layers = layers, n_latent = latents)
    vae_adata_red.train(early_stopping=True, max_epochs = 400)
    print(colored('Training complete.', 'cyan'))
    # plot learning curve
    train_elbo = vae_adata_red.history['elbo_train'][1:]
    test_elbo= vae_adata_red.history['elbo_validation']
    ax = train_elbo.plot()
    test_elbo.plot(ax=ax)
    plt.rcParams['figure.figsize'] = 5,5
    plt.show(block=False)
    adata_latent = vae_adata_red.get_latent_representation()
    adata.obsm['X_scVI'] = adata_latent
    sc.pp.neighbors(adata, use_rep='X_scVI')
    sc.tl.umap(adata)
    res = simpledialog.askfloat('Leiden resolution parameter', 'Type desired leiden resolution')
    sc.tl.leiden(adata, resolution = res)
    #check in query
    sns.set(rc={'figure.figsize':(10,10)})
    sns.set_style('whitegrid')
    sc.pl.umap(adata, color=['leiden'], size=40, legend_loc='on data')
    plt.show(block=False)
    answer = messagebox.askyesno('Keep resolution?', message= 'Do you wish to keep the resolution parameter?')
    while answer == False:
        res = simpledialog.askfloat('Leiden resolution parameter', 'Type desired leiden resolution')
        sc.tl.leiden(adata, resolution = res)
        sc.pl.umap(adata, color=['leiden'], size=40, legend_loc="on data")
        plt.show(block=False)
        answer = messagebox.askyesno('Keep resolution?', message= 'Do you wish to keep the resolution parameter?')
    print(colored('-------------------------------------------------------------------------\n', 'white'), colored('Workflow complete! You can now continue with plotting and quantification.', 'green'))
    # declare object as global for use in plotting and writing functions, also make copy to recall after subetting
    adata.obs['predictions'] = adata.obs['leiden']
    global adata_to_sub
    global bdata
    adata_to_sub = adata.copy()
    bdata = adata.copy()
    bdata.layers['scvi_norm'] = vae_adata_red.get_normalized_expression(library_size=10e4)

def ask_hypers():
    hyper_window = Toplevel()
    hyper_window.geometry('250x200')
    hyper_window.title('Choose hyperparameters for Neural Network')
    Message(hyper_window, text = 'Number of neurons per layer (default = 128):', width=350).pack()
    hyper1 = Entry(hyper_window, textvariable = num_neurons)
    hyper1.pack()
    Message(hyper_window, text = 'Number of layers (default = 1):', width=350).pack()
    hyper2 = Entry(hyper_window, textvariable = num_layers)
    hyper2.pack()
    Message(hyper_window, text = 'Number of latent variables (default = 10)', width=350).pack()
    hyper3 = Entry(hyper_window, textvariable = num_latents)
    hyper3.pack()
    # make apply button
    submit_var = IntVar()
    submit = Button(hyper_window, text = "Train model", command = lambda:[submit_var.set(1), hyper_window.destroy()])
    submit.pack()
    submit.wait_variable(submit_var)

def ask_threshs():
    thresh_window = Toplevel()
    thresh_window.geometry('200x150')
    thresh_window.title('Choose thresholds for QC variables')
    Message(thresh_window, text = 'Min counts:', width=350).pack()
    thresh1 = Entry(thresh_window, textvariable = count_thresh)
    thresh1.pack()
    Message(thresh_window, text = 'Min unique features:', width=350).pack()
    thresh2 = Entry(thresh_window, textvariable = feat_thresh)
    thresh2.pack()
    Message(thresh_window, text = 'Max %MT:', width=350).pack()
    thresh3 = Entry(thresh_window, textvariable = mt_thresh)
    thresh3.pack()
    # make apply button
    submit = Button(thresh_window, text = "Apply thresholds", command = lambda:[thresh_window.destroy()])
    submit.pack()

# ask for .tsv file holding condition identifier
def load_condition_csv():
	cond = filedialog.askopenfilename(title = 'Load condition.tsv')
	path_condition.set(value = cond)
	print(path_condition.get())
	global condition
	path = path_condition.get()
	condition = pd.read_csv(path, sep=' ')

# ask for .tsv file holding sample identifier
def load_ident_csv():
    ident = filedialog.askopenfilename(title = 'Load ident.tsv')
    path_ident.set(value = ident)
    print(path_ident.get())
    global idents
    path = path_ident.get()
    idents = pd.read_csv(path, sep=' ')

# ask for .tsv file holding study identifier
def load_study_csv():
	study = filedialog.askopenfilename(title = 'Load study.tsv')
	path_study.set(value = study)
	print(path_study.get())
	global study_ident
	path = path_study.get()
	study_ident = pd.read_csv(path, sep=' ')

# checkbuttons
def cbutt_freq():
    check1 = check1.get()
    cbutt1

# make dge functions
def dge_betw_celltypes():
    print(colored('Differentially expressed genes between cell types will be computed...', 'yellow'))
    dge_data = bdata.copy()
    # get chosen subset
    subset = subset_choice5.get()
    if subset == 'All cells':
        dge_data_sub = dge_data.copy()
    elif subset == 'CD4 T cells':
        dge_data_sub = dge_data[dge_data.obs.predictions.isin(['Tfh', 'Th17', 'Th2/Th22', 'Th1', 'Tregs', 'CCR5high Th17.1', 'CD4 TEMRA'])].copy()
    elif subset == 'CD8 T cells':
        dge_data_sub = dge_data[dge_data.obs.predictions.isin(['MAIT', 'gdT Vd2+', 'gdT Vd2-', 'CD8 CM', 'CD8 EM HLA-DRA+', 'CD8 EM CD160+', 'CD8 TRM ITGA1+', 'CD8 TRM ITGA1-', 'CD8 CTL'])].copy()
    elif subset == 'Myeloid cells':
        dge_data_sub = dge_data[dge_data.obs.predictions.isin(['Monocytes', 'BAM MRC1+', 'BAM EMP3+', 'MG CX3CR1+', 'MG CCL2+', 'MG TREM2hi', 'mDCs CD1c+', 'mDCs AXL+SIGLEC6+', 'mDCs CLEC9A+'])].copy()
    elif subset == 'NK cells':
        dge_data_sub = dge_data[dge_data.obs.predictions.isin(['NK bright', 'NK dim', 'TR-NK', 'ILC'])].copy()
    elif subset == 'B cells':
        dge_data_sub = dge_data[dge_data.obs.predictions.isin(['B cells', ' B activated', 'B IL4R+', 'B atypical', 'Plasmacells'])].copy()
    # dge
    dge_model = vae_adata_red.differential_expression(dge_data_sub, groupby = 'predictions', batch_correction=True)
    # make marker list
    markers = {}
    cats = pd.Series(dge_data_sub.obs['predictions'].values, dtype="category").cat.categories
    for i, c in enumerate(cats):
        cid = "{} vs Rest".format(c)
        cell_type_df = dge_model.loc[dge_model.comparison == cid]
        cell_type_df = cell_type_df[cell_type_df.lfc_mean > 0]
        cell_type_df = cell_type_df[cell_type_df["bayes_factor"] > 0]
        cell_type_df = cell_type_df[cell_type_df["non_zeros_proportion1"] > 0.1]
        markers[c] = cell_type_df.index.tolist()[:10]
    sc.pl.matrixplot(dge_data_sub, markers, cmap="inferno", groupby="predictions", standard_scale='var', dendrogram=True, layer = 'scvi_norm')
    make_stats_window(x=dge_model)
    reset_dge3()

def dge_betw_conditions():
    dge_data = bdata.copy()
    # get chosen subset
    subset = subset_choice4.get()
    print(colored('Differentially expressed genes between conditions for ', 'yellow') + subset + colored(' will be computed...'))
    # make subset
    if subset == 'All cells':
        dge_data_sub = dge_data.copy()
    elif subset == 'CD4 T cells':
        dge_data_sub = dge_data[dge_data.obs.predictions.isin(['Tfh', 'Th17', 'Th2/Th22', 'Th1', 'Tregs', 'CCR5high Th17.1', 'CD4 TEMRA'])].copy()
    elif subset == 'CD8 T cells':
        dge_data_sub = dge_data[dge_data.obs.predictions.isin(['MAIT', 'gdT Vd2+', 'gdT Vd2-', 'CD8 CM', 'CD8 EM HLA-DRA+', 'CD8 EM CD160+', 'CD8 TRM ITGA1+', 'CD8 TRM ITGA1-', 'CD8 CTL'])].copy()
    elif subset == 'Myeloid cells':
        dge_data_sub = dge_data[dge_data.obs.predictions.isin(['Monocytes', 'BAM MRC1+', 'BAM EMP3+', 'MG CX3CR1+', 'MG CCL2+', 'MG TREM2hi', 'mDCs CD1c+', 'mDCs AXL+SIGLEC6+', 'mDCs CLEC9A+'])].copy()
    elif subset == 'NK cells':
        dge_data_sub = dge_data[dge_data.obs.predictions.isin(['NK bright', 'NK dim', 'TR-NK', 'ILC'])].copy()
    elif subset == 'B cells':
        dge_data_sub = dge_data[dge_data.obs.predictions.isin(['B cells', ' B activated', 'B IL4R+', 'B atypical', 'Plasmacells'])].copy()
    else:
        dge_data_sub = dge_data[dge_data.obs.predictions.isin([subset])].copy()
    # dge
    dge_model = vae_adata_red.differential_expression(dge_data_sub, groupby = 'condition', batch_correction=True)
    # make marker list
    markers = {}
    cats = pd.Series(dge_data_sub.obs['condition'].values, dtype="category").cat.categories
    for i, c in enumerate(cats):
        cid = "{} vs Rest".format(c)
        cell_type_df = dge_model.loc[dge_model.comparison == cid]
        cell_type_df = cell_type_df[cell_type_df.lfc_mean > 0]
        cell_type_df = cell_type_df[cell_type_df["bayes_factor"] > 0]
        cell_type_df = cell_type_df[cell_type_df["non_zeros_proportion1"] > 0.1]
        markers[c] = cell_type_df.index.tolist()[:10]
    sc.pl.matrixplot(dge_data_sub, markers, cmap="inferno", groupby="condition", standard_scale='var', dendrogram=True, layer = 'scvi_norm')
    make_stats_window(x=dge_model)
    reset_dge2()

# destroy subwindow
def close(x):
    x.destroy()

# make subwindow for
def subwindow3():
    quant_win2 = Toplevel()

    # get logo
    global image2
    if getattr(sys, 'frozen', False):
        image2 = PIL.Image.open(os.path.join(sys._MEIPASS, "files/Plot_options.png"))
    else:
        image2 = PIL.Image.open("files/Plot_options.png")

    #image = image.resize((143*8,430))
    image2 =  PIL.ImageTk.PhotoImage(image2)
    Label(quant_win2, image = image2).pack(side=TOP)

    quant_win2.title('Which plot do you want to create?')
    # checkbuttons
    Button(quant_win2 , text = 'Stacked barchart', width = 30, height = 1, command = lambda:[plot_freq1()]).pack(side=LEFT)
    Button(quant_win2 , text = 'Boxplot by celltype', width = 54, height = 1, command = lambda:[plot_freq2()]).pack(side=LEFT)
    Button(quant_win2, text = 'Boxplot by celltype and condition', width = 56, height = 1, command = lambda:[plot_freq3()]).pack(side=LEFT)
    heatmap_options = ['Default lineage/subset heatmap', 'Custom heatmap']
    OptionMenu(quant_win2, heatmap_choice, *heatmap_options, command = lambda _:[plot_heatmap()]).pack(side=LEFT)

def subwindow4():
    dge_win = Toplevel()
    # image for dge
    try:
        global image3
        if getattr(sys, 'frozen', False):
            image3 = PIL.Image.open(os.path.join(sys._MEIPASS, "files/dge_options.png"))
        else:
            image3 = PIL.Image.open("files/dge_options.png")

        #image = image.resize((143*8,430))
        image3 =  PIL.ImageTk.PhotoImage(image3)
        Label(dge_win, image = image3).pack(side=TOP)
        # Button for "Between cells"
        OptionMenu(dge_win, subset_choice5, *options5, command = lambda _:[dge_betw_celltypes()]).pack(side=LEFT)
        # Dropdown for "Between condition"
        OptionMenu(dge_win, subset_choice4, *options4, command = lambda _:[dge_betw_conditions()]).pack(side=LEFT)
    except:
        print(colored('These subsets are only defined for the mapping workflow. Choose "All cells"'))

def keep_all():
    subset_choice2.set('Do not subset')
    subset_choice.set('Do not subset')

def basic_workflow_window():
    bflow = Toplevel()
    # get logo
    global image4
    # get logo
    if getattr(sys, 'frozen', False):
        image4 = PIL.Image.open(os.path.join(sys._MEIPASS, "files/acsf_logo_basic_workflow.png"))
    else:
        image4 = PIL.Image.open("files/acsf_logo_basic_workflow.png")
    image4 = image4.resize((1300,200))
    image4 =  PIL.ImageTk.PhotoImage(image4)
    wl = Label(bflow, image = image4)
    wl.pack()
    # make buttons
    Button(bflow, text = 'I. a. Choose 10X data', command = lambda:[get_path()], width = 19, height = 1).pack(side=LEFT)
    Button(bflow, text = 'II. a. Choose condition.tsv', command = lambda:[load_condition_csv()], width = 19, height = 1).pack(side=LEFT)
    Button(bflow, text = 'II. b. Choose study.tsv', command = lambda:[load_study_csv()], width = 19, height = 1).pack(side=LEFT)
    Button(bflow, text = 'II. c. Choose idents.tsv', command = lambda:[load_ident_csv()], width = 19, height = 1).pack(side=LEFT)
    Button(bflow, text= 'III. Run basic workflow', command = lambda:[load_user_adata2()], width = 19, height = 1).pack(side=LEFT)
    Button(bflow, text= 'IV. a. Plot UMAP', command = lambda :[get_genes_to_plot(), keep_all(), plot_umap(), recall_full_data()], width = 19, height = 1).pack(side=LEFT)
    Button(bflow, text = 'IV. b. Differential expression', command = subwindow4).pack(side=LEFT)
    Button(bflow, text = 'VI. Rename leiden clusters', command = lambda: [get_new_names(), map_new_names(), keep_all()]).pack(side=LEFT)

def get_new_names():
    clusters = bdata.obs['predictions'].value_counts().index.tolist()
    global string_vars
    string_vars = {}
    for i in range(0, len(clusters)):
        locals()['clust' + str(i)] = StringVar()
        string_vars[i] = locals()['clust' + str(i)]
    rename_window = Toplevel()
    for i in range(0, len(clusters)):
        text_output = 'New name for clust' + str(i)
        Message(rename_window, text = 'Type new name for cluster ' + str(i), width=350).pack()
        locals()['rename_', str(i)] = Entry(rename_window, textvariable = string_vars[i])
        locals()['rename_', str(i)].pack()
    submit = Button(rename_window, text = "Apply thresholds", command = lambda:[submit_var.set(1), rename_window.destroy()])
    submit_var = IntVar()
    submit.pack()
    submit.wait_variable(submit_var)

def map_new_names():
    # map to predictions
    print(bdata)
    print('start_loop')
    clusters = bdata.obs['predictions'].value_counts().index.tolist()
    for i in range(0, len(clusters)):
        print(bdata)
        print(colored('Changing name of Cluster...', 'yellow'), colored(i, 'white'))
        bdata.obs['predictions'] = np.where(bdata.obs['predictions'] == str(i), string_vars[i].get(), bdata.obs['predictions'])
    print(bdata.obs['predictions'].value_counts())

# make var for dir
path_var = StringVar(root)
path_condition = StringVar(root)
path_study = StringVar(root)
path_ident = StringVar(root)
subset_choice = StringVar(root)
subset_choice.set('IV. b. Plot quantifications')
subset_choice3 = StringVar(root)
subset_choice3.set('IV. c. Differential expression')
subset_choice5 = StringVar(root)
subset_choice5.set('Between cell types')
subset_choice4 = StringVar(root)
subset_choice4.set('Between conditions')
data_to_write = StringVar(root)
data_to_write.set('IV. d. Write data')
subset_choice2 = StringVar(root)
subset_choice2.set('IV. a. Plot UMAP')
heatmap_choice = StringVar(root)
heatmap_choice.set('Heatmap for label markers')
count_thresh = IntVar(root)
count_thresh.set(1500)
feat_thresh = IntVar(root)
feat_thresh.set(500)
mt_thresh = IntVar(root)
mt_thresh.set(4)
num_neurons = IntVar(root)
num_neurons.set(128)
num_layers = IntVar(root)
num_layers.set(1)
num_latents = IntVar(root)
num_latents.set(10)
learnrate = IntVar(root)


# make button to start loading scvi data
butt1 = Button(root, text = 'I. a. Choose 10X data', command = lambda:[get_path()], width = 19, height = 1)
butt2 = Button(root, text = 'II. a. Choose condition.tsv', command = lambda:[load_condition_csv()], width = 19, height = 1)
butt3 = Button(root, text = 'II. b. Choose study.tsv', command = lambda:[load_study_csv()], width = 19, height = 1)
butt4 = Button(root, text = 'II. c. Choose idents.tsv', command = lambda:[load_ident_csv()], width = 19, height = 1)
butt5 = Button(root, text= 'III. Run mapping', command = lambda:[load_user_adata()], width = 19, height = 1)
options3 = ['Do not subset', 'CD4 T cells', 'CD8 T cells', 'Myeloid cells', 'NK cells', 'B cells']
butt6 = OptionMenu(root, subset_choice2,  *options3, command = lambda _:[get_genes_to_plot(), plot_umap(), reset_plot_umap(), recall_full_data()])
options = ['Do not subset', 'CD4 T cells', 'CD8 T cells', 'Myeloid cells', 'NK cells', 'B cells']
butt7 = OptionMenu(root, subset_choice, *options, command = lambda _:[subset_fun(), reset_quant(), subwindow3(), recall_full_data()])
options2 = ['Annotations', 'UMAP coordinates', '.h5ad']
butt8 = OptionMenu(root, data_to_write, *options2, command = lambda _:[write_data(), reset_write()])
butt9 = Button(root, text = 'IV. c. Differential expression', command = subwindow4)
global options4
options4 = ['All cells', 'CD4 T cells', 'CD8 T cells', 'Myeloid cells', 'NK cells', 'B cells', ' B activated', 'B IL4R+', 'B atypical', 'Plasmacells', 'pDCs', 'Monocytes', 'BAM MRC1+', 'BAM EMP3+', 'MG CX3CR1+', 'MG CCL2+', 'MG TREM2hi', 'mDCs CD1c+', 'mDCs AXL+SIGLEC6+', 'mDCs CLEC9A+', 'NK bright', 'NK dim', 'TR-NK', 'ILC', 'MAIT', 'gdT Vd2+', 'gdT Vd2-', 'CD8 CM', 'CD8 EM HLA-DRA+', 'CD8 EM CD160+', 'CD8 TRM ITGA1+', 'CD8 TRM ITGA1-', 'CD8 CTL', 'Tfh', 'Th17', 'Th2/Th22', 'Th1', 'Tregs', 'CCR5high Th17.1', 'CD4 TEMRA']
global options5
options5 = ['All cells', 'CD4 T cells', 'CD8 T cells', 'Myeloid cells', 'NK cells', 'B cells']
butt10 = Button(root, text= 'Basic workflow', command = basic_workflow_window, width = 19, height = 1)

# make var to store genes to plot
gene_var = StringVar(root)
order_var = StringVar(root)
hue_var = StringVar(root)
stats_adata_var = StringVar(root)


def reset_quant():
    subset_choice.set('IV. b. Plot quantifications')

def reset_dge():
    subset_choice3.set('IV. c. Differential expression')

def reset_dge2():
    subset_choice4.set('Between conditions')

def reset_dge3():
    subset_choice5.set('Between cell types')

def reset_write():
    data_to_write.set('IV. d. Write data')

def reset_plot_umap():
    subset_choice2.set('IV. a. Plot UMAP')

# add labels to widget
wl.pack(side=TOP)
butt1.pack(side=LEFT)
butt2.pack(side=LEFT)
butt3.pack(side=LEFT)
butt4.pack(side=LEFT)
butt5.pack(side=LEFT)
butt6.pack(side=LEFT)
butt7.pack(side=LEFT)
butt9.pack(side=LEFT)
butt8.pack(side=LEFT)
butt10.pack(side=LEFT)

# run app
root.mainloop()
