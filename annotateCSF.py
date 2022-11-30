#!/usr/bin/env python

#tk imports
from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
from tkinter import simpledialog
from tkinter import ttk
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
root.title('annotate CSF')

# welcome statement
def welcome_statement():
    print(colored("\nWelcome to aCSF! \nYou can use aCSF for label transfer and subsequent plotting and analysis of your CSF scRNAseq data. Please refer to the reference manual for more detailed information regarding the usage of aCSF.\nFor any issues with the software please raise an issue on our github repo: github.com/schwab-lab/annotateCSF\n\nPatrick Ostkamp, 2022 \n", 'yellow'), colored("\nThis tool was built using scanpy, anndata, scVI, and scANVI workflows. When using, please consider citing the original publications as well: \n\nWolf et al. (2018), Scanpy: large-scale single-cell gene expression data analysis, Genome Biology. \n\nVirshup et al. (2021) anndata: Annotated data, bioRxiv. \n\nGayoso, Lopez, Xing, et al (2022), A Python library for probabilistic analysis of single-cell omics data, Nature Biotechnology. \n\nLopez (2018), Deep generative modeling for single-cell transcriptomics, Nature Methods. \n\nXu, Lopez, Mehlman, et al (2021), Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models, Molecular Systems Biology.\n-----------------------------------------------------------------------------------------------------------------------\n", "white" ))
welcome_statement()

# get logo
if getattr(sys, 'frozen', False):
    image = PIL.Image.open(os.path.join(sys._MEIPASS, "files/acsf_logo.png"))
else:
    image = PIL.Image.open("files/acsf_logo.png")
# add logo to root
image = image.resize((1480,560))
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

# get order of conditions for box-/jitter plot from user
def get_condition_order():
    cond_order_window = Toplevel()
    # get list of conditions
    conditions = list(adata_copy.obs['condition'].cat.categories.values)
    print(conditions)
    cond_order = []
    for i in range(0,len(conditions)):
        label_text = "Choose condition no. " + str((i+1))
        Label(cond_order_window, text = label_text).pack()
        wait_var = StringVar()
        def activate():
            wait_var.set(1)
            cond_box.config(state='disabled')
        chosen_cond = StringVar()
        cond_box = ttk.Combobox(cond_order_window, values = conditions, textvariable = chosen_cond)
        cond_box.pack()
        cond_butt = Button(cond_order_window, text = "Confirm", command = lambda:[activate()])
        cond_butt.pack()
        cond_butt.wait_variable(wait_var)
        cond_order.append(chosen_cond.get())
    cond_order_window.destroy()
    global condition_order
    condition_order = cond_order

# get hues to color conditions accordingly in box-/jitter plot
def get_hues():
    cols = {}
    col_prompts = list(['Choose first color', 'Choose second color', 'Choose third color', 'Choose fourth color', 'Choose fifth color', 'Choose sixth color', 'Choose seventh color', 'Choose eighth color', 'Choose nineth color', 'Choose tenth'])
    print(len(adata_copy.obs.condition.value_counts().index.tolist()))
    for i in range(0, len(adata_copy.obs.condition.value_counts().index.tolist())):
        col_choice = askcolor(title=col_prompts[i])
        cols[i] = str(col_choice[1])
    hue = list(cols.values())
    hue_var.set(value=hue)
    global new_hues
    new_hues = hue_var

# prepare data
def prep_data():
    # prepare widget
    prep = Toplevel()
    Label(prep, text = "Choose directory with folders holding 10X data").grid(row=1,column=1)
    # get user dir with 10X folders
    def get_dir():
        user_dir = filedialog.askdirectory(title="Select directory containing folders with 10X data. \nThis should only contain folders holding matrix.mtx, barcodes.tsv, and genes.tsv or features.tsv")
        dir.set(user_dir)
    dir = StringVar()
    dir_button = Button(prep, text = "Select directory", textvariable = dir, command = lambda:[get_dir()])
    dir_button.grid(row=2,column=1)
    dir_button.wait_variable(dir)
    # list of folders in dir
    folders = os.listdir(dir.get())
    print("Chosen directory is:" + dir.get())
    no_of_folders = len(folders)
    # define get_sample_info
    sample_list =[0] * no_of_folders
    study_list =[0] * no_of_folders
    condition_list =[0] * no_of_folders

    print(no_of_folders, "folders were provided")
    def get_sample_info():
        for i in range(0,no_of_folders):
            # make variables to fill
            Label(prep, text = "Define metadata for sample: " + folders[i]).grid(row=4,column=1)
            sample_list[i] = StringVar()
            study_list[i] = StringVar()
            condition_list[i] = StringVar()
            Label(prep, text = "Sample identifier (orig.ident):").grid(row=5,column=1)
            Label(prep, text = "Study identifier:").grid(row=5,column=3)
            Label(prep, text = "Condition identifier:").grid(row=5,column=5)
            # make entry fields
            sample_name = Entry(prep, textvariable = sample_list[i])
            sample_name.grid(row=5, column=2)
            study_name = Entry(prep, textvariable = study_list[i])
            study_name.grid(row=5, column=4)
            condition_name = Entry(prep, textvariable = condition_list[i])
            condition_name.grid(row=5, column=6)
            # wait for confirm
            wait_var_loop = IntVar()
            def activate():
                wait_var_loop.set(1)
            def enable_finish():
                if i == (no_of_folders-1):
                    finish.config(state='normal')
            submit = Button(prep, text = "Confirm", command = lambda:[activate(), enable_finish()])
            submit.grid(row=5,column=7)
            submit.wait_variable(wait_var_loop)
            print("done with " + folders[i] + "!")

    # confirm folder
    submit = Button(prep, text = "Confirm", command = lambda:[get_sample_info()])
    submit.grid(row = 2, column = 2)

    # wait for get_sample_info to be finished
    wait_sample_info = IntVar()
    def activate_finish():
        wait_sample_info.set(1)
    finish = Button(prep, text = "Start data preparation", command = lambda:[activate_finish()])
    finish.grid(row=6, column = 1)
    finish.config(state="disabled")
    finish.wait_variable(wait_sample_info)

    # load actual data
    datasets = []
    print("Starting to load data...")
    for i in range(0,no_of_folders):
        # get gene matrix
        file = os.listdir((dir.get() + "/" + folders[i]))[2]
        path_to_mat = dir.get() + "/" + folders[i] + "/" + file
        print(path_to_mat)
        adata = sc.read(path_to_mat).T
        #check whether features or genes.tsv is provided
        file = os.listdir((dir.get() + "/" + folders[i]))[1]
        genes_path = dir.get() + "/" + folders[i]  + "/" + file
        # read genes
        genes = pd.read_csv(genes_path, header=None, sep='\t')
        adata.var['gene_symbols'] = genes[1].values
        adata.var_names = adata.var['gene_symbols']
        adata.var_names_make_unique(join="-")
        # load barcodes
        file = os.listdir((dir.get() + "/" + folders[i]))[0]
        barcodes_path = dir.get() + "/" + folders[i] + "/" + file
        cells = pd.read_csv(barcodes_path, header=None, sep='\t')
        adata.obs['barcode'] = cells[0].values
        adata.obs_names = cells[0]
        adata.obs_names_make_unique(join="-")
        # add user input as metadata
        adata.obs['orig.ident'] = sample_list[i].get()
        adata.obs['study'] = study_list[i].get()
        adata.obs['condition'] = condition_list[i].get()
        datasets.append(adata)
    print(datasets)
    prep.destroy()

    # concatenate data
    print("Concatenating data...")
    adata = datasets[0]
    for i in range(1,len(datasets)):
        adata = adata.concatenate(datasets[i])
    print(adata)

    ## write data
    messagebox.showinfo(title="Writing data...", message = "Choose a directory to store the 10X formatted data")
    # write 10x data
    write_dir = filedialog.askdirectory(title = "Where to store the 10X formatted data?")
    from scipy.sparse import csr_matrix
    import scipy.io as sio
    final_mat = csr_matrix(adata.X.T)
    path = write_dir + "/" + "matrix.mtx"
    sio.mmwrite(path, final_mat)
    path = write_dir + "/" + "barcodes.tsv"
    barcodes = pd.Series(adata.obs_names.values)
    barcodes.to_csv(path, index=False, header=False)
    path = write_dir + "/" + "features.tsv"
    features = pd.Series(adata.var_names.values)
    features.to_csv(path, index=False, header=False)

    # write metadata
    messagebox.showinfo(title="Writing data...", message = "Choose a directory to store the metadata")
    write_dir = filedialog.askdirectory(title = "Where to store the metadata?")
    path = write_dir + "/" + "idents.tsv"
    adata.obs['orig.ident'].to_csv(path, index=False)
    path = write_dir + "/" + "study.tsv"
    adata.obs['study'].to_csv(path, index=False)
    path = write_dir + "/" + "condition.tsv"
    adata.obs['condition'].to_csv(path, index=False)

# Compute enrichment scores
def enrich_score():
    # interface
    enrich_window = Toplevel()
    Label(enrich_window, text = "Compute enrichment score for predefined gene expression signatures\nPrepare your signature as single .tsv files containing only gene symbols (also no header)\n").pack()
    signature_dir = StringVar()
    def ask_signature_dir():
         dir = filedialog.askopenfilename(title="Choose gene expression signature .tsv file")
         signature_dir.set(dir)
    chosen_score_name = StringVar()
    Label(enrich_window, text = "Enter a name for the signature:").pack()
    Entry(enrich_window, textvariable = chosen_score_name).pack()
    dir_button = Button(enrich_window, text = "Choose .tsv file", command = lambda:[ask_signature_dir()])
    dir_button.pack()
    wait_var = StringVar()
    confirm_button = Button(enrich_window, text = "Confirm", command = lambda:[wait_var.set(1)])
    confirm_button.pack()
    confirm_button.wait_variable(wait_var)
    enrich_window.destroy()
    chosen_score = pd.read_csv(signature_dir.get(), names = ['score'])
    # how many genes were not found?
    number_in_sig = len(chosen_score['score'])
    number_found_in_adata = len(set(chosen_score['score']).intersection(set(adata_copy.raw.var_names)))
    number_call = str(number_found_in_adata) + " genes of " + str(number_in_sig) + " supplied were found."
    messagebox.showinfo(message = number_call)
    # compute enrichment
    sc.tl.score_genes(adata_copy, chosen_score['score'], score_name = chosen_score_name.get(), use_raw = True)
    # make available for umap etc.
    if "chosen_scores" not in globals():
        global chosen_scores
        chosen_scores = list([chosen_score_name.get()])
    else:
        chosen_scores.append(chosen_score_name.get())

# plot default heatmap
def plot_heatmap_def():
    plot1 = {' B activated | B IL4R+ | B atypical | Plasmacells': ['CD19', 'CD79A', 'TNFRSF13B'], 'pDCs': ['IL3RA', 'IRF8', 'LAMP5'], 'Monocytes | BAM MRC1+ | BAM EMP3+ | MG CX3CR1+ | MG CCL2+ | MG TREM2hi': ['CD14', 'CD68', 'MS4A7'], 'mDCs CD1c+ | mDCs AXL+SIGLEC6+ | mDCs CLEC9A+': ['CD1C', 'AFF3', 'HLA-DRB5'], 'NK bright | NK dim | TR-NK | ILC': ['NCAM1', 'GNLY', 'XCL1'], 'MAIT | gdT Vd2+ | gdT Vd2- | CD8 CM | CD8 EM HLA-DRA+ | CD8 EM CD160+ | CD8 TRM ITGA1+ | CD8 TRM ITGA1- | CD8 CTL': ['CD8B', 'CD8A', 'GZMH'], 'Tfh | Th17 | Th2/Th22 | Th1 | Tregs | CCR5high Th17.1 | CD4 TEMRA': ['CD4', 'TNFRSF25', 'AQP3']}
    plot2 = {' B activated': ['CD69', 'CXCR4', 'BACH2'], 'B IL4R+': ['IL4R', 'FCER2', 'IGHM'], 'B atypical': ['FCRL5', 'CD1C', 'GPR34'], 'Plasmacells': ['IGHG1', 'SDC1', 'CD38'], 'pDCs': ['IL3RA', 'IRF7', 'IRF8'], 'Monocytes': ['VCAN', 'S100A12', 'CCR2'], 'BAM MRC1+': ['MRC1', 'KCNAB1', 'MARCO'], 'BAM EMP3+': ['EMP3', 'CYP27A1', 'TIMD4'], 'MG CX3CR1+': ['CX3CR1', 'TMEM119', 'P2RY12'], 'MG CCL2+': ['SPP1', 'IRAK2', 'ITGAX'], 'MG TREM2hi': ['TREM2', 'APOC1', 'GPNMB'], 'mDCs CD1c+': ['CD1C', 'CD207', 'CD1E'], 'mDCs AXL+SIGLEC6+': ['SIGLEC6', 'CLIC3', 'CD5D'], 'mDCs CLEC9A+': ['CLEC9A', 'PPP1R14A', 'TMEM14A'], 'NK bright': ['NCAM1', 'KLRC3', 'SPINK2'], 'NK dim': ['EOMES', 'TTC38', 'FGFBP2'], 'TR-NK': ['IKZF3', 'KRT81', 'ITGAE'], 'ILC': ['KIT', 'TNFRSF4', 'IFNGR2'], 'MAIT': ['TRAV1-2', 'NCR3', 'KLRB1'], 'gdT Vd2+': ['TRDC', 'TRGV9', 'ZBTB16'], 'gdT Vd2-': ['TRDV2', 'LSR', 'RTKN2'], 'CD8 CM': ['CCR7', 'NELL2', 'MAL'], 'CD8 EM HLA-DRA+': ['HLA-DRA', 'MSC', 'HLA-DRB5'], 'CD8 EM CD160+': ['CD160', 'FCRL6', 'FGR'], 'CD8 TRM ITGA1+': ['ITGA1', 'ZNF683', 'ITGAE'], 'CD8 TRM ITGA1-': ['NR4A2', 'CEMIP2', 'CXCR4'], 'CD8 CTL': ['GNLY', 'FXYD2', 'GZMB'], 'Th17.1,CD4': ['TEMRA', 'PASK', 'CXCR5'], 'Th17': ['CCR6', 'USP10', 'RORC'], 'Th2/Th22': ['SLC40A1', 'SOCS1', 'CCR4'], 'Th1': ['TBX21', 'TRBC1', 'CXCR3'], 'Tregs': ['FOXP3', 'RTKN2', 'IKZF2'], 'CCR5high Th17.1': ['CXCR6', 'CD69', 'RGS1'], 'CD4 TEMRA': ['GZMH', 'PDCD1', 'CX3CR1']}

    subsets_in_data = set(adata_sub2.obs.predictions.value_counts().index.values)   # gibt eine Liste der vorhandenen subsets im dataset
    print(adata_sub2.raw.var)
    genes_in_data_set = set(adata_sub2.raw.var.index.values) # extrahiert die Liste von vorhandenen Genen aus dem dataset

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

# plot custom heatmap
def plot_heatmap_cust():
    # genes
    cust_heatmap_window = Toplevel()
    available_genes = list(adata_sub2.raw.var.index.values)
    chosen_genes = []
    sorted(available_genes)
    label_text = "Choose genes to plot"
    Label(cust_heatmap_window, text = label_text).grid(row = 0, column = 1)
    for i in range(0,999):
        chosen_gene = StringVar()
        wait_var_loop = StringVar()
        choose_genes = ttk.Combobox(cust_heatmap_window, values = available_genes, state = "choose", textvariable = chosen_gene)
        choose_genes.grid(row = 1, column = 1)
        def activate():
            wait_var_loop.set(2)
        confirm_button = Button(cust_heatmap_window, text = "Confirm choice", command= lambda:[activate()])
        confirm_button.grid(row = 1, column = 2)
        confirm_button.wait_variable(wait_var_loop)
        chosen_genes.append(chosen_gene.get())
        label_text = "The following genes will be plotted:"
        Label(cust_heatmap_window, text = label_text).grid(row=2, column = 1)
        Label(cust_heatmap_window, text =  str(chosen_genes)).grid(row=3, column = 1)
        choose_genes.configure(state="disabled")
        wait_plot = IntVar()
        wait_plot.set(1)
        more_button = Button(cust_heatmap_window, text = "Add another gene", command = lambda:[wait_plot.set(1)])
        plot_button = Button(cust_heatmap_window, text = "Finish and Plot", command = lambda:[wait_plot.set(2)])
        more_button.grid(row = 1, column = 3)
        plot_button.grid(row = 1, column = 4)
        plot_button.wait_variable(wait_plot)
        plot_button.configure(state='disabled')
        more_button.configure(state='disabled')
        confirm_button.configure(state='disabled')
        if wait_plot.get() == 2:
            cust_heatmap_window.destroy()
            break
    #submit_button.wait_variable(wait_var)
    sc.pl.matrixplot(adata_sub2, var_names= chosen_genes, groupby='predictions', standard_scale='var', dendrogram=True, cmap='inferno')

# decide custom or default heatmap
def plot_heatmap():
    if (heatmap_choice.get() == 'Custom heatmap'):
        plot_heatmap_cust()
    else:
        plot_heatmap_def()

# plot umap
def plot_umap():
    ## Receive data
    # get chosen subset
    subset = subset_choice2.get()

    ## User choices
    # make window
    cb_genes = Toplevel()
    Label(cb_genes, text = "I. Choose data to plot.").grid(column = 1, row = 1)

    # colors, size, etc
    # palettes
    pal_cont_var = StringVar()
    pal_cont_var.set("viridis")
    pal_cat_var =  StringVar()
    pal_cat_var.set("glasbey_light")
    # define glasby_light palette
    glasb = ['#d60000', '#018700', '#b500ff', '#05acc6', '#97ff00', '#ffa52f', '#ff8ec8', '#79525e', '#00fdcf', '#afa5ff', '#93ac83', '#9a6900', '#366962', '#d3008c', '#fdf490', '#c86e66', '#9ee2ff', '#00c846', '#a877ac', '#b8ba01', '#f4bfb1', '#ff28fd', '#f2cdff', '#009e7c', '#ff6200', '#56642a', '#953f1f', '#90318e', '#ff3464', '#a0e491', '#8c9ab1', '#829026', '#ae083f', '#77c6ba', '#bc9157']
    glasbey_light = sns.color_palette(glasb)
    # available palettes
    pals_cont = ['viridis', 'inferno', 'magma', 'turbo', 'mako', 'plasma', 'cividis', 'rocket']
    pals_cat = ['glasbey_light', 'turbo', 'Pastel1', 'Pastel2', 'Paired', 'Accent', 'Dark2', 'Set1', 'Set2', 'Set3', 'tab10', 'tab20', 'tab20b', 'tab20c']
    pal_cont_label = Label(cb_genes, text = "III. Choose continuous color palette")
    pal_cont_label.grid(column = 4, row = 1)
    choose_pal_cont = ttk.Combobox(cb_genes, values = pals_cont, state = "viridis", textvariable = pal_cont_var)
    choose_pal_cont.grid(column = 4, row = 2)
    pal_cat_label = Label(cb_genes, text = "IV. Choose categorical color palette")
    pal_cat_label.grid(column = 4, row = 4)
    choose_pal_cat = ttk.Combobox(cb_genes, values = pals_cat, state = "glasbey_light", textvariable = pal_cat_var)
    choose_pal_cat.grid(column =4, row = 5)
    # point size and legend size
    pt_size = StringVar()
    pt_size.set(50)
    pt_size_label = Label(cb_genes, text = "V. Choose point size")
    pt_size_label.grid(column = 6, row = 1)
    pt_size_entry = Entry(cb_genes, textvariable = pt_size)
    pt_size_entry.grid(column = 6, row = 2)
    lg_size = StringVar()
    lg_size.set(10)
    lg_size_label = Label(cb_genes, text = "VI. Choose legend fontsize")
    lg_size_label.grid(column = 6, row = 4)
    lg_size_entry = Entry(cb_genes, textvariable = lg_size)
    lg_size_entry.grid(column = 6, row = 5)

    # define values for combobox
    metadata = ['predictions', 'lineage', 'predictions_l2', 'condition', 'orig.ident', 'study'] # all available metadata
    all_genes = list(adata_copy.raw.var.index.values) # available genes in dataset
    all_genes.sort() # alphabetical
    available_info = metadata + all_genes
    # add gene scores if computed before...
    if "chosen_scores" in globals():
        available_info = chosen_scores + available_info

    # genes
    chosen_genes = StringVar()
    choose_genes = ttk.Combobox(cb_genes, values = available_info, state = "choose", textvariable = chosen_genes)
    # define genes_to_plot variable upon selection
    def create_genes_to_plot(self):
        global genes_to_plot
        chosen = chosen_genes.get()
        genes_to_plot = list()
        genes_to_plot.append(chosen)
        print("Chosen to plot:", genes_to_plot)
    # add combobox to widget
    choose_genes.bind("<<ComboboxSelected>>", create_genes_to_plot)
    choose_genes.grid(column = 1, row = 2)

    # make button to execute plot. Waits for activation.
    activate_plot = IntVar()
    plot_button = Button(cb_genes, text = "PLOT!", command = lambda:[activate_plot.set(1)])
    plot_button.grid(column = 6, row = 6)

    # Add more info to plot. Requires genes_to_plot to be defined
    def add_info():
        try:
            # combobox to add further info
            chosen_genes_add = StringVar()
            choose_genes_add = ttk.Combobox(cb_genes, values = available_info, state = "predictions", textvariable = chosen_genes_add)
            choose_genes_add.grid(column = 1, row = 6)
            # confirm adding by pressing button
            confirm_dummy = DoubleVar()
            def set_confirm_dummy():
                confirm_dummy.set(random.uniform(0,9999999))
            confirm_add_button = Button(cb_genes, text = "Confirm", command = set_confirm_dummy)
            confirm_add_button.grid(column = 1, row = 7)
            # wait until button is pressed
            confirm_add_button.wait_variable(confirm_dummy) # waits for pressing confirm_add_button. Button generates random float, therefore variable will always get changed upon button press to stop wait_variable command
            # after button press, get chosen gene
            choose_genes_add.config(state="disabled") # disable until add_more_info is pressed again
            chosen = chosen_genes_add.get()
            genes_to_plot.append(chosen)
        except:
            print("Select a medadata column or gene from the dropdown menu first")

    add_info_label = Label(cb_genes, text="II. Choose further data and hit confirm")
    add_info_label.grid(column = 1, row = 4)
    add_info_button = Button(cb_genes, text = "Add more info", command = lambda:[add_info()])
    add_info_button.grid(column = 1, row = 5)

    # stop waiting...
    plot_button.wait_variable(activate_plot)

    ## Execute plot
    if subset == 'Do not subset':
        try:
            if pal_cat_var.get() == "glasbey_light":
                sc.pl.umap(adata_copy, color = genes_to_plot, size=int(pt_size.get()), palette=glasbey_light, cmap=pal_cont_var.get(), legend_loc="on data", legend_fontsize = int(lg_size.get()))
            else:
                sc.pl.umap(adata_copy, color = genes_to_plot, size=int(pt_size.get()), palette=pal_cat_var.get(), cmap=pal_cont_var.get(), legend_loc="on data", legend_fontsize = int(lg_size.get()))
        except:
            print("Some error occured...")

    # Define chosen subset
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

    # If subset was not plotted before, perform scanpy or scVI workflow to generate subset UMAP
    if (subset != 'Do not subset') and ('adata_sub' not in globals()):
        print(colored('You chose to subset your data for the UMAP plot. To provide better resolution, a basic scanpy workflow will be run for your chosen subset.', 'yellow'))
        global adata_sub
        adata_sub = adata_copy.copy()
        print("this is adata_sub:", adata_sub)
        # this needs to be simplified: Lineage will be defined in mapping() and should be usable here.
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
        # subset adata to only contain one lineage
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
        print(pt_size.get(), pal_cat_var.get(), pal_cont_var.get(), lg_size.get())
        sc.pl.umap(adata_sub, color = genes_to_plot, size=int(pt_size.get()), palette=pal_cat_var.get(), cmap=pal_cont_var.get(), legend_loc="on data", legend_fontsize = int(lg_size.get()))
        try:
            print("This should be the subset plot...")
            if pal_cat_var.get() == "glasbey_light":
                sc.pl.umap(adata_sub, color = genes_to_plot, size=int(pt_size.get()), palette=glasbey_light, cmap=pal_cont_var.get(), legend_loc="on data", legend_fontsize = int(lg_size.get()))
            else:
                sc.pl.umap(adata_sub, color = genes_to_plot, size=int(pt_size.get()), palette=pal_cat_var.get(), cmap=pal_cont_var.get(), legend_loc="on data", legend_fontsize = int(lg_size.get()))
        except:
            #print('You entered ' + genes_to_plot + '. This is either not a valid gene symbol or the gene was not found in the dataset')
            print("Some error occured...")
    # If workflow has been run for subset already re-use that data
    else:
        if (subset != 'Do not subset') and ('adata_sub' in globals()) and (subset in adata_sub.obs.mnc_lineage.values):
            try:
                if pal_cat_var.get() == "glasbey_light":
                    sc.pl.umap(adata_sub, color = genes_to_plot, size=int(pt_size.get()), palette=glasbey_light, cmap=pal_cont_var.get(), legend_loc="on data", legend_fontsize = int(lg_size.get()))
                else:
                    sc.pl.umap(adata_sub, color = genes_to_plot, size=int(pt_size.get()), palette=pal_cat_var.get(), cmap=pal_cont_var.get(), legend_loc="on data", legend_fontsize = int(lg_size.get()))
            except:
                #print('You entered ' + genes_to_plot + '. This is either not a valid gene symbol or the gene was not found in the dataset')
                print("Some error occured...")

    if 'adata_sub' in globals():
        print('adata_sub is in globals...')
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
        adata_sub = adata_copy.copy()
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
        try:
            if pal_cat_var.get() == "glasbey_light":
                sc.pl.umap(adata_sub, color = genes_to_plot, size=int(pt_size.get()), palette=glasbey_light, cmap=pal_cont_var.get(), legend_loc="on data", legend_fontsize = int(lg_size.get()))
            else:
                sc.pl.umap(adata_sub, color = genes_to_plot, size=int(pt_size.get()), palette=pal_cat_var.get(), cmap=pal_cont_var.get(), legend_loc="on data", legend_fontsize = int(lg_size.get()))
        except:
            print("Gene or metadata not available...")

# subset to mnc
def keep_mncs():
    global adata_sub2
    adata_sub2 = adata_copy.copy()

# subset function
def subset_fun():
    global adata_sub2
    subset = subset_choice.get()
    if subset == 'CD4 T cells':
        adata_sub2 = adata_copy[adata_copy.obs.predictions.isin(['Tfh', 'Th17', 'Th2/Th22', 'Th1', 'Tregs', 'CCR5high Th17.1', 'CD4 TEMRA'])]
    elif subset == 'CD8 T cells':
        adata_sub2 = adata_copy[adata_copy.obs.predictions.isin(['CD8 TRM ITGA1+', 'CD8 TRM ITGA1-', 'CD8 EM HLA-DRA+', 'CD8 EM CD160+', 'CD8 CM', 'CD8 CTL', 'gdT Vd2+', 'gdT Vd2-', 'MAIT'])]
    elif subset == 'Myeloid cells':
        adata_sub2 = adata_copy[adata_copy.obs.predictions.isin(['MG CX3CR1+', 'Monocytes', 'MG CCL2+', 'MG TREM2hi', 'BAM MRC1+', 'BAM EMP3+', 'mDCs CD1c+', 'mDCs AXL+SIGLEC6+', 'mDCs CLEC9A+'])]
    elif subset == 'NK cells':
        adata_sub2 = adata_copy[adata_copy.obs.predictions.isin(['NK bright', 'NK dim', 'TR-NK', 'ILC'])]
    elif subset == 'B cells':
        adata_sub2 = adata_copy[adata_copy.obs.predictions.isin(['B IL4R+', ' B activated', 'Plasmacells', 'B atypical'])]
    elif subset == 'Do not subset':
        adata_sub2 = adata_copy.copy()

# recall full anndata after subsetting
def recall_full_data():
    global adata_copy
    adata_copy = adata_copy.copy()

# stacked barplot
def stacked_bar_plot():
    anot_df = adata_sub2.obs['predictions'].value_counts()
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

# box-/jitterplot by celltype
def boxplot_by_celltype():
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

# Stats functions
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
    celltypes = quant_df.celltype.value_counts().index
    for i in range(0,len(celltypes)):
      df = quant_df[quant_df.celltype==celltypes[i]]
      form = "percentage ~ C(condition, Treatment(reference = '" + stats_adata_group + "'))"
      model = ols(form, data = df)
      fii = model.fit()
      lm_res[i] = fii.summary2().tables[1]
      lm_res[i]['celltype'] = celltypes[i]
      conds = df.condition.value_counts().index.categories.to_list()
      conds.remove(stats_adata_group)
      lm_res[i] = lm_res[i].iloc[1:,:]
      comparison = pd.DataFrame(list([stats_adata_group])*len(conds)) + pd.DataFrame(list([' vs '])*len(conds)) + pd.DataFrame(conds)
      comparison[0].values
      lm_res[i]['comparison'] = comparison[0].values
    lm_res2 = pd.concat(lm_res)
    make_stats_window(x=lm_res2)

# define boxplot by celltype and condition
def boxplot_by_condition():

    # make df from adata
    anot_df = adata_sub2.obs.copy()

    # define boxplot scaffold
    def box_jitter(which_df, cell_order, pal, cond_order):
        sns.set(rc={'figure.figsize':(30,7)})
        sns.set_style(rc={'grid.color':'0.9', 'axes.facecolor':'white', 'axes.edgecolor':'grey'})
        ax = sns.boxplot(x = 'celltype', y='percentage', data = which_df, hue='condition', palette = len(anot_df.condition.value_counts())*["white"], order=cell_order, hue_order = cond_order, dodge=True, fliersize = 0, showmeans=True, meanprops={'marker':'+', 'markersize':'6', 'markeredgewidth':'1.5', 'markeredgecolor':'#777777'})
        ax = sns.stripplot(x= 'celltype', y='percentage', data = which_df, alpha=.5, hue = 'condition', dodge=True, order = cell_order, palette = pal, hue_order = cond_order)
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode='anchor')
        ax.tick_params(bottom=True, left=True)
        plt.ylabel('% of sample')
        plt.show()

    # define customization
    def customize():
        # new condition order
        get_condition_order()
        new_order = condition_order
        # new colors
        get_hues()
        hues = new_hues.get()
        hues = list(eval(hues))
        colors = hues
        return new_order, colors

    # declare global variable for stats function
    global quant_df

    # define set of all possible subsets
    all_cells = [' B activated', 'B IL4R+', 'B atypical', 'Plasmacells', 'pDCs', 'Monocytes', 'BAM MRC1+', 'BAM EMP3+', 'MG CX3CR1+', 'MG CCL2+', 'MG TREM2hi', 'mDCs CD1c+', 'mDCs AXL+SIGLEC6+', 'mDCs CLEC9A+', 'NK bright', 'NK dim', 'TR-NK', 'ILC', 'MAIT', 'gdT Vd2+', 'gdT Vd2-', 'CD8 CM', 'CD8 EM HLA-DRA+', 'CD8 EM CD160+', 'CD8 TRM ITGA1+', 'CD8 TRM ITGA1-', 'CD8 CTL', 'Tfh', 'Th17', 'Th2/Th22', 'Th1', 'Tregs', 'CCR5high Th17.1', 'CD4 TEMRA']
    all_cells_lin = ['CD4 T cells', 'CD8 T cells', 'B cells', 'NK cells', 'Myeloid cells']

    # if data was not subsetted, show lineage plot first, then all subsets
    if len(anot_df.lineage.value_counts()) > 1:
        # aggregate data to percentages
        lin_df = anot_df.copy()
        lin_df['count'] = 1
        variables = pd.DataFrame({'orig.ident':anot_df['orig.ident'], 'condition':anot_df['condition']}).drop_duplicates()
        lin_df = lin_df.groupby(['orig.ident', 'lineage', 'condition'],sort=False).agg({'count':'sum'}).groupby(level=0).apply(lambda x: 100*x/x.sum()).groupby(level = ['orig.ident', 'lineage']).max().reset_index()
        lin_df = lin_df.merge(variables, left_on = 'orig.ident', right_on = 'orig.ident', how='left').drop_duplicates()
        lin_df['count'].sum()
        lin_df = lin_df.set_axis(['orig.ident', 'celltype', 'percentage', 'condition'], axis=1)
        # make set of lineages found in the dataset and find intersection with all possible lineages
        available_cells = set(lin_df.celltype.value_counts().index)
        print(all_cells_lin)
        print(available_cells)
        cell_order = list(set(all_cells_lin).intersection(set(available_cells)))
        print(cell_order)
        # default lineage plot
        pal = 'turbo'
        default_order = lin_df.condition.value_counts().index.values
        box_jitter(which_df = lin_df, cell_order = cell_order, pal = pal, cond_order = default_order)
        # store quantification results for lineages?
        answer = messagebox.askyesno(title="Store results", message = "Do you want to store quantification results?")
        if answer == True:
            path = filedialog.askdirectory(title = 'Where to store quantification results?') + "/lineage_quantification.csv"
            lin_df.to_csv(path)
        # store statistics?
        answer = messagebox.askyesno(title="Store results", message = "Generate statistics with linear model?")
        if answer == True:
            quant_df = lin_df
            get_stats_adata()
            make_stats()
        # customize plot?
        answer = messagebox.askyesno(title="Customize plot?", message = "Do you want to customize the plot?")
        if answer == True:
            new_order, colors = customize()
            box_jitter(which_df = lin_df, cell_order = cell_order, pal = colors, cond_order=new_order)

    # if data was subsetted, only plot subset
    # aggregate to percentages
    anot_df['count'] = 1
    variables = pd.DataFrame({'orig.ident':anot_df['orig.ident'], 'condition':anot_df['condition']}).drop_duplicates()
    anot_df = anot_df.groupby(['orig.ident', 'predictions', 'condition'],sort=False).agg({'count':'sum'}).groupby(level=0).apply(lambda x: 100*x/x.sum()).groupby(level = ['orig.ident', 'predictions']).max().reset_index()
    anot_df = anot_df.merge(variables, left_on = 'orig.ident', right_on = 'orig.ident', how='left').drop_duplicates()
    anot_df['count'].sum()
    anot_df = anot_df.set_axis(['orig.ident', 'celltype', 'percentage', 'condition'], axis=1)
    # find intersection of found subsets and all possible subsets
    available_cells = set(anot_df.celltype.value_counts().index)
    cell_order = list(set(all_cells).intersection(set(available_cells)))
    default_order = anot_df.condition.value_counts().index.values
    # default plot with all subsets
    box_jitter(which_df = anot_df, cell_order = cell_order, pal = 'turbo', cond_order= default_order)
    # store quantification results for lineages?
    answer = messagebox.askyesno(title="Store results?", message = "Do you want to store quantification results?")
    if answer == True:
        path = filedialog.askdirectory(title = 'Where to store quantification results?') +"/" + str(anot_df.lineage.value_counts().index) + "_quantification.csv"
        anot_df.to_csv(path)
    # store statistics?
    answer = messagebox.askyesno(title="Store results", message = "Generate statistics with linear model?")
    if answer == True:
        quant_df = anot_df
        get_stats_adata()
        make_stats()
    # customize plot?
    answer = messagebox.askyesno(title="Customize plot?", message = "Do you want to customize the plot?")
    if answer == True:
        new_order, colors = customize()
        print(colors)
        box_jitter(which_df = anot_df, cell_order = cell_order, pal = colors, cond_order=new_order)

# write data to disk: h5ad, annotations, umap coords
def write_data():
    data_type = data_to_write.get()
    get_path()
    if data_type == '.h5ad':
        print(colored('Writing to .h5ad...', 'red'))
        #adata_copy.var.index.name = "gene symbols"
        from scipy.sparse import csr_matrix
        adata_copy.X = csr_matrix(adata_copy.X)
        path = path_var.get()
        path = path + "/labelled_query.h5ad"
        print('Writing data to:' + path)
        adata_copy.write_h5ad(path)
        print(colored('...done!', 'blue'))
    elif data_type == 'Annotations':
        print(colored('Writing annotations to .tsv...', 'red'))
        csv_out = pd.DataFrame({'mapped_labels':adata_copy.obs.predictions, 'orig.ident':adata_copy.obs['orig.ident'], 'study':adata_copy.obs['study'], 'condition':adata_copy.obs['condition']})
        path = path_var.get()
        path = path + "/mapped_labels.tsv"
        print('Writing data to:' + path)
        csv_out.to_csv(path, sep='\t')
        print(colored('...done!', 'cyan'))
    elif data_type == 'UMAP coordinates':
        print(colored('Writing umap coords to .tsv...', 'red'))
        umap_out = pd.DataFrame(adata_copy.obsm['X_umap'])
        umap_out = umap_out.rename({0: 'UMAP_1', 1: 'UMAP_2'}, axis=1)
        path = path_var.get()
        path = path + "/umap_coords.tsv"
        umap_out.to_csv(path, sep='\t')
        print(colored('...done!', 'cyan'))

# main scANVI workflow
def map_to_query():

    print(colored('Loading 10X data...', 'red'))
    # load mtx file
    matrix_path = path_var.get() + '/matrix.mtx'
    global adata
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
    ## ! This is a quick solution and may need adjustment depending on
    genes = pd.read_csv(genes_path, header=None, sep='\t')
    if genes[0].str.contains('ENSG000').any() == False:
        try:
            adata.var['gene_symbols'] = genes[0].values
        except:
            adata.var['gene_symbols'] = genes[1].values
    if genes[0].str.contains('ENSG000').any() == True:
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
        adata.obs['orig.ident'] = idents.values
    else:
    	adata.obs['orig.ident'] = "sample"

    if 'condition' in globals():
    	adata.obs['condition'] = condition.values
    else:
    	adata.obs['condition'] = "condition"

    if 'study_ident' in globals():
    	adata.obs['study'] = study_ident.values
    else:
    	adata.obs['study'] = "study"

    # find doublets with scrublet
    answer = messagebox.askyesno('Doublet discimination', 'Perform doublet discrimination with scrublet?')
    if answer == True:
        print(colored('Finding doublets...', 'red'))
        batch_list = []
        # make vector of sample names
        samples = adata.obs["orig.ident"].unique()
        for x in samples:
            print('starting', x)
            batch = adata[adata.obs["orig.ident"].isin([x])]
            scrub = scr.Scrublet(batch.layers["counts"])
            doublet_scores, predicted_doublets = scrub.scrub_doublets()
            print('doublets found')
            batch.obs['doublets'] = predicted_doublets
            batch.obs['doublet_score'] = doublet_scores
            print('added to object')
            batch_list.append(batch)
            print(colored('Batch', 'cyan'), colored(x, 'white'), colored(' done!', 'cyan'))
        # merge datasets
        adata = batch_list[0]
        for x in range(1, (len(batch_list))):
            adata = adata.concatenate(batch_list[x])
        plt.hist(adata.obs.doublet_score, bins = 50)
        plt.show(block = False)
        doublet_thresh = simpledialog.askfloat(title = "Doublet threshold", prompt = "Type desired doublet threshold (max value)")
        plt.close()
        adata = adata[adata.obs.doublet_score < doublet_tresh]
        print(colored('Doublets found!', 'cyan'))

    # normalization for HVG calculation
    print(colored('Preprocessing and QC...', 'red'))
    sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed = True)
    sc.pp.log1p(adata)
    adata.raw = adata

    # QC vars
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['rp'] = adata.var_names.str.startswith('RPS', 'RPL')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'rp'], percent_top=None, inplace=True)
    # no of cells before QC
    print(adata)

    # to wait until subsetting has been done...
    global wait_subset
    wait_subset = IntVar()

    # ask how to qc
    qc_window = Toplevel()
    qc_label = Label(qc_window, text = "How to perform QC?\n Sample by sample (based on idents.tsv) is recommended for data from several batches.")
    qc_label.pack()
    qc_decide_sbs = Button(qc_window, text = "Sample by sample", command = lambda:[ask_threshs_sbs(), subset_sbs(), qc_window.destroy()])
    qc_decide_sbs.pack()
    qc_decide_ato = Button(qc_window, text = "All samples at once", command = lambda:[ask_threshs(), subset_ato(), qc_window.destroy()])
    qc_decide_ato.pack()

    # wait until wait_subset is set by subset_sbs or subset_ato
    root.wait_variable(wait_subset)

    # show number of removed cells
    print(adata)

    # load Ostkamp et al. labeled reference
    print(colored('Loading reference data...', 'red'))
    if getattr(sys, 'frozen', False):
        ref = sc.read_h5ad(os.path.join(sys._MEIPASS, "files/ref2.h5ad"))
    else:
        ref = sc.read_h5ad("files/ref2.h5ad")
    print(colored('Reference loaded!', 'cyan'))

    # assign batch for discrimination after merging
    ref.obs['batch'] = '0'
    adata.obs['batch'] = '1'

    # rename CD8 TRM labels
    ref.obs['labels'] = np.where(ref.obs['labels']=='CD8 TRM1', 'CD8 TRM ITGA1+', ref.obs['labels'])
    ref.obs['labels'] = np.where(ref.obs['labels']=='CD8 TRM2', 'CD8 TRM ITGA1-', ref.obs['labels'])
    ref.obs['labels'] = np.where(ref.obs['labels']=='CD8 proliferating', 'Proliferating', ref.obs['labels'])


    # preprocess reference
    sc.pp.normalize_total(ref, target_sum=1e4, exclude_highly_expressed = True)
    sc.pp.log1p(ref)
    ref.raw = ref

    print(ref)
    print(ref.var.index)
    print(adata.var.index)

    # combine, find hvgs, and split
    print(colored('Concatenating datasets and finding variable genes...', 'red'))
    combined = ref.concatenate(adata)
    print(combined)
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

    # set seed to random integer
    random.seed()
    scvi.settings.seed = random.randint(0,99999999)

    # make VAE to pre-train reference
    print(colored('Preparing data and commencing training on reference dataset', 'red'))
    scvi.data.setup_anndata(ref_red, layer="counts", batch_key="orig.ident", categorical_covariate_keys = ['tissue', 'study'], continuous_covariate_keys = ['pct_counts_mt'], labels_key='labels')
    ref_model = scvi.model.SCVI(ref_red, n_hidden = neurons[0], n_layers = 2, n_latent = latents[0], dropout_rate=.2)
    ref_model.train(early_stopping=True, max_epochs = 20)
    print(colored('Training complete.', 'cyan'))
    # obtain latent and compute umap for reference
    ref_latent = ref_model.get_latent_representation()
    ref_red.obsm['X_scVI'] = ref_latent
    sc.pp.neighbors(ref_red, use_rep='X_scVI')
    sc.tl.umap(ref_red)
    # train scANVI VAE for ref
    print(colored('Supervising trained model using labels...', 'red'))
    vae_scanvi_model = scvi.model.SCANVI.from_scvi_model(ref_model, unlabeled_category='Unknown', adata=ref_red)
    vae_scanvi_model.train(max_epochs = 50, early_stopping = True)
    print(colored('...done!', 'cyan'))
    ref_red.obsm['X_scANVI'] = vae_scanvi_model.get_latent_representation()
    sc.pp.neighbors(ref_red, use_rep='X_scANVI')
    sc.tl.umap(ref_red)
    # set labels
    adata_red.obs['labels'] = 'Unknown'
    ref_red.obs['predictions'] = vae_scanvi_model.predict()

    # map labels
    global vae_adata_red
    vae_adata_red = scvi.model.SCANVI.load_query_data(adata_red, vae_scanvi_model)
    print(colored('Starting mapping process... This may take a while. Using GPU via cuda may accelerate the process.', 'red'))
    vae_adata_red.train(max_epochs = 400, early_stopping = True, check_val_every_n_epoch = 1, plan_kwargs=dict(weight_decay=0.0), early_stopping_min_delta =0.01, early_stopping_patience = 15)
    adata_red.obsm['X_scANVI'] = vae_adata_red.get_latent_representation()
    adata_red.obs['predictions'] = vae_adata_red.predict()
    sc.pp.neighbors(adata_red, use_rep = 'X_scANVI')
    sc.tl.umap(adata_red)
    print(colored('-------------------------------------------------------------------------\n', 'white'), colored('Mapping complete! You can now continue with plotting and quantification.', 'green'))

    # define clustering levels
    sub_to_lin = dict({'Th2/Th22':'CD4 T cells', 'Tfh':'CD4 T cells', 'Th1':'CD4 T cells', 'MG CX3CR1+':'Myeloid cells', 'CD8 EM HLA-DRA+':'CD8 T cells', 'CCR5high Th17.1':'CD4 T cells', 'CD4 TEMRA':'CD4 T cells', 'CD8 TRM ITGA1-':'CD8 T cells', 'CD8 CM':'CD8 T cells', 'BAM EMP3+':'Myeloid cells', 'Th17':'CD4 T cells', 'Tregs':'CD4 T cells', 'BAM MRC1+':'Myeloid cells', 'CD8 EM CD160+':'CD8 T cells', 'CD8 TRM ITGA1+':'CD8 T cells', 'mDCs CD1c+': 'Myeloid cells', 'MG CCL2+':'Myeloid cells', 'MG TREM2hi':'Myeloid cells', 'Monocytes':'Myeloid cells', 'CD8 CTL': 'CD8 T cells', 'pDCs':'pDCs', 'NK bright':'NK cells', 'gdT Vd2+':'CD8 T cells',' B activated':'B cells', 'NK dim':'NK cells', 'B atypical':'B cells', 'gdT Vd2-':'CD8 T cells', 'TR-NK':'NK cells', 'MAIT':'CD8 T cells','ILC':'NK cells','Plasmacells':'B cells', 'mDCs CLEC9A+':'Myeloid cells', 'B IL4R+':'B cells', 'mDCs AXL+SIGLEC6+':'Myeloid cells', 'CD8 prolif':'CD8 T cells'})
    adata_red.obs['lineage'] = adata_red.obs['predictions'].replace(sub_to_lin)
    # l2
    sub_to_l2 = dict({'Th2/Th22':'CD4 helper T', 'Tfh':'CD4 helper T', 'Th1':'CD4 helper T', 'MG CX3CR1+':'Microglia', 'CD8 EM HLA-DRA+':'CD8 effector T', 'CCR5high Th17.1':'CD4 helper T', 'CD4 TEMRA':'CD4 TEMRA', 'CD8 TRM ITGA1-':'CD8 effector T', 'CD8 CM':'CD8 CM', 'BAM EMP3+':'BAM', 'Th17':'CD4 helper T', 'Tregs':'CD4 Tregs', 'BAM MRC1+':'BAM', 'CD8 EM CD160+':'CD8 effector T', 'CD8 TRM ITGA1+':'CD8 effector T', 'mDCs CD1c+': 'mDCs', 'MG CCL2+':'Microglia', 'MG TREM2hi':'Microglia', 'Monocytes':'Monocytes', 'CD8 CTL': 'CD8 effector T', 'pDCs':'pDCs', 'NK bright':'NK cells', 'gdT Vd2+':'gamma-delta T',' B activated':'B cells', 'NK dim':'NK cells', 'B atypical':'B cells', 'gdT Vd2-':'gamma-delta T', 'TR-NK':'TR-NK', 'MAIT':'MAIT','ILC':'ILC','Plasmacells':'Plasmacells', 'mDCs CLEC9A+':'mDCs', 'B IL4R+':'B cells', 'mDCs AXL+SIGLEC6+':'mDCs', 'CD8 prolif':'Proliferating T'})
    adata_red.obs['predictions_l2'] = adata_red.obs['predictions'].replace(sub_to_l2)

    # declare object as global for use in plotting and writing functions, also make copy to recall after subsetting
    global adata_copy
    # define adata_copy as new working data, remove adata_red
    adata_copy = adata_red.copy()
    del(adata_red)
    adata_copy.layers['scvi_norm'] = vae_adata_red.get_normalized_expression(library_size=10e4)

def scVI_workflow():

    print(colored('Loading 10X data...', 'red'))
    # load mtx file
    matrix_path = path_var.get() + '/matrix.mtx'
    global adata
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
    try:
        adata.var['gene_symbols'] = genes[1].values
    except:
        adata.var['gene_symbols'] = genes[0].values
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
        adata.obs['orig.ident'] = idents.values
    else:
    	adata.obs['orig.ident'] = "sample"

    if 'condition' in globals():
    	adata.obs['condition'] = condition.values
    else:
    	adata.obs['condition'] = "condition"

    if 'study_ident' in globals():
    	adata.obs['study'] = study_ident.values
    else:
    	adata.obs['study'] = "study"

    # find doublets with scrublet
    answer = messagebox.askyesno('Doublet discimination', 'Perform doublet discrimination with scrublet?')
    if answer == True:
        print(colored('Finding doublets...', 'red'))
        batch_list = []
        # make vector of sample names
        samples = adata.obs["orig.ident"].unique()
        for x in samples:
            print('starting', x)
            batch = adata[adata.obs["orig.ident"].isin([x])]
            scrub = scr.Scrublet(batch.layers["counts"])
            doublet_scores, predicted_doublets = scrub.scrub_doublets()
            print('doublets found')
            batch.obs['doublets'] = predicted_doublets
            batch.obs['doublet_score'] = doublet_scores
            print('added to object')
            batch_list.append(batch)
            print(colored('Batch', 'cyan'), colored(x, 'white'), colored(' done!', 'cyan'))
        # merge datasets
        adata = batch_list[0]
        for x in range(1, (len(batch_list))):
            adata = adata.concatenate(batch_list[x])
        plt.hist(adata.obs.doublet_score, bins = 50)
        plt.show(block = False)
        doublet_thresh = simpledialog.askfloat(title = "Doublet threshold", prompt = "Type desired doublet threshold (max value)")
        plt.close()
        adata = adata[adata.obs.doublet_score < 0.2]
        print(colored('Doublets found!', 'cyan'))

    # normalization for HVG calculation
    print(colored('Preprocessing and QC...', 'red'))
    sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed = True)
    sc.pp.log1p(adata)
    adata.raw = adata

    # QC vars
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    adata.var['rp'] = adata.var_names.str.startswith('RPS', 'RPL')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'rp'], percent_top=None, inplace=True)
    # no of cells before QC
    print(adata)

    # to wait until subsetting has been done...
    global wait_subset
    wait_subset = IntVar()

    # ask how to qc
    qc_window = Toplevel()
    qc_label = Label(qc_window, text = "How to perform QC?\n Sample by sample (based on idents.tsv) is recommended for data from several batches.")
    qc_label.pack()
    qc_decide_sbs = Button(qc_window, text = "Sample by sample", command = lambda:[ask_threshs_sbs(), subset_sbs(), qc_window.destroy()])
    qc_decide_sbs.pack()
    qc_decide_ato = Button(qc_window, text = "All samples at once", command = lambda:[ask_threshs(), subset_ato(), qc_window.destroy()])
    qc_decide_ato.pack()

    # wait until wait_subset is set by subset_sbs or subset_ato
    root.wait_variable(wait_subset)

    # preprocess adataerece
    sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed = True)
    sc.pp.log1p(adata)
    adata.raw = adata
    # combine, find hvgs, and split
    print(colored('Finding variable genes...', 'red'))
    # ask user number of hvgs
    n_hvgs = simpledialog.askinteger(title='Highly variable genes', prompt='Enter desired number of highly variable genes')
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvgs, subset=True, flavor = 'seurat', batch_key='orig.ident')
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
    global vae_adata
    vae_adata = scvi.model.SCVI(adata, n_hidden = neurons, n_layers = layers, n_latent = latents)
    vae_adata.train(early_stopping=True, max_epochs = 400)
    print(colored('Training complete.', 'cyan'))
    # plot learning curve
    train_elbo = vae_adata.history['elbo_train'][1:]
    test_elbo= vae_adata.history['elbo_validation']
    ax = train_elbo.plot()
    test_elbo.plot(ax=ax)
    plt.rcParams['figure.figsize'] = 5,5
    plt.show(block=False)
    # get latent representation and compute umap
    adata_latent = vae_adata.get_latent_representation()
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
    # declare working objects
    global adata_copy
    # make copy
    adata_copy = adata.copy()
    adata_copy.layers['scvi_norm'] = vae_adata.get_normalized_expression(library_size=10e4)
    # remove adata
    del(adata)

# get hyperparameters for neural network
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

# get user qc thresholds for all samples at once
def ask_threshs():
    # show qc plot
    global adata
    plt.rcParams['figure.figsize'] = 8,8
    fig, axes = plt.subplots(2,2)
    sns.violinplot(data=adata.obs, x="condition", y="total_counts", ax=axes[0,0], color="white")
    sns.stripplot(data=adata.obs, x="condition", y="total_counts", ax=axes[0,0], alpha = 0.15, size = 3)
    sns.violinplot(data=adata.obs, x="condition", y="n_genes_by_counts", ax=axes[0,1], color="white")
    sns.stripplot(data=adata.obs, x="condition", y="n_genes_by_counts", ax=axes[0,1], alpha=0.15, size = 3)
    sns.violinplot(data=adata.obs, x="condition", y="pct_counts_mt", ax=axes[1,0], color="white")
    sns.stripplot(data=adata.obs, x="condition", y="pct_counts_mt", ax=axes[1,0], alpha=0.15, size = 3)
    sns.scatterplot(data=adata.obs, x="total_counts", y="n_genes_by_counts", ax=axes[1,1], hue =  "pct_counts_mt", palette = 'viridis', size=3)
    plt.tight_layout()
    plt.show(block=False)
    def plot_close():
        plt.close()
    # ask threshold
    thresh_window = Toplevel()
    thresh_window.geometry('250x150+0+0')
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
    print("abcd")
    wait_var = StringVar()
    def activate():
        wait_var.set(1)
    # make apply button
    submit = Button(thresh_window, text = "Apply thresholds", command = lambda:[activate(), thresh_window.destroy()])
    submit.pack()
    submit.wait_variable(wait_var)
    plot_close()
    print("Thresholds defined!")

# get user qc thresholds one by one
def ask_threshs_sbs():
    global adata
    global count_threshs
    count_threshs = []
    global feat_threshs
    feat_threshs = []
    global mt_threshs
    mt_threshs = []
    for i in range(0,len(adata.obs['orig.ident'].value_counts().index.values)):
        # make qc plot
        plt.rcParams['figure.figsize'] = 8,8
        fig, axes = plt.subplots(2,2)
        sns.violinplot(data=adata[adata.obs['orig.ident'] == adata.obs['orig.ident'].value_counts().index.values[i]].obs, x="condition", y="total_counts", ax=axes[0,0], color="white")
        sns.stripplot(data=adata[adata.obs['orig.ident'] == adata.obs['orig.ident'].value_counts().index.values[i]].obs, x="condition", y="total_counts", ax=axes[0,0], alpha = 0.25, size = 3)
        sns.violinplot(data=adata[adata.obs['orig.ident'] == adata.obs['orig.ident'].value_counts().index.values[i]].obs, x="condition", y="n_genes_by_counts", ax=axes[0,1], color="white")
        sns.stripplot(data=adata[adata.obs['orig.ident'] == adata.obs['orig.ident'].value_counts().index.values[i]].obs, x="condition", y="n_genes_by_counts", ax=axes[0,1], alpha=0.25, size = 3)
        sns.violinplot(data=adata[adata.obs['orig.ident'] == adata.obs['orig.ident'].value_counts().index.values[i]].obs, x="condition", y="pct_counts_mt", ax=axes[1,0], color="white")
        sns.stripplot(data=adata[adata.obs['orig.ident'] == adata.obs['orig.ident'].value_counts().index.values[i]].obs, x="condition", y="pct_counts_mt", ax=axes[1,0], alpha=0.25, size = 3)
        sns.scatterplot(data=adata[adata.obs['orig.ident'] == adata.obs['orig.ident'].value_counts().index.values[i]].obs, x="total_counts", y="n_genes_by_counts", ax=axes[1,1], hue =  "pct_counts_mt", palette = 'viridis', size=3)
        plt.tight_layout()
        plt.show(block=False)
        def plot_close():
            plt.close()
        # ask for thresholds
        threshs_window = Toplevel()
        identity = "QC for batch: " + adata.obs['orig.ident'].value_counts().index.values[i]
        threshs_label = Label(threshs_window, text = identity)
        threshs_label.pack()
        threshs_window.geometry('350x200+0+0')
        threshs_window.title('Choose qc thresholds for batch')
        Message(threshs_window, text = 'Min counts:', width=350).pack()
        threshs1 = Entry(threshs_window, textvariable = count_thresh)
        threshs1.pack()
        Message(threshs_window, text = 'Min unique features:', width=350).pack()
        threshs2 = Entry(threshs_window, textvariable = feat_thresh)
        threshs2.pack()
        Message(threshs_window, text = 'Max %MT:', width=350).pack()
        threshs3 = Entry(threshs_window, textvariable = mt_thresh)
        threshs3.pack()
        def get_threshs():
            count_threshs.append(count_thresh.get())
            feat_threshs.append(feat_thresh.get())
            mt_threshs.append(mt_thresh.get())
        # make apply button
        wait_var = StringVar()
        def activate():
            wait_var.set(1)
        # make apply button
        submit = Button(threshs_window, text = "Apply thresholds", command = lambda:[get_threshs(), activate(), threshs_window.destroy()])
        submit.pack()
        submit.wait_variable(wait_var)
        plot_close()

# subset according to user qc thresholds
# at once
def subset_ato():
    global adata
    print("Subsetting data according to thresholds...")
    adata = adata[(adata.obs.total_counts > count_thresh.get()) & (adata.obs.n_genes_by_counts > feat_thresh.get()) & (adata.obs.pct_counts_mt < mt_thresh.get())]
    print(colored('...done!', 'cyan'))
    wait_subset.set(1)

# sample by sample
def subset_sbs():
    global adata
    global count_threshs
    global feat_threshs
    global mt_threshs
    global wait_subset
    print(count_threshs)
    batch_list = []
    batches = adata.obs['orig.ident'].value_counts().index.values
    print(batches)
    for i in range(0,len(adata.obs['orig.ident'].value_counts().values)):
        adata_batch = adata[adata.obs['orig.ident']==batches[i]]
        adata_batch = adata_batch[(adata_batch.obs.total_counts > count_threshs[i]) & (adata_batch.obs.n_genes_by_counts > feat_threshs[i]) & (adata_batch.obs.pct_counts_mt < mt_threshs[i])]
        batch_list.append(adata_batch)
    if len(batch_list) == 1:
        adata_concat = batch_list[0]
    else:
        for j in range(1, (len(batch_list))):
            adata_concat = batch_list[0].concatenate(batch_list[j])
    print("Dataset after QC contains:", adata_concat.n_obs)
    adata = adata_concat.copy()
    wait_subset.set(1)
    del(adata_concat)

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
    dge_data = adata_copy.copy()
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
    try:
        dge_model = vae_adata_red.differential_expression(dge_data_sub, groupby = 'predictions')
    except:
        dge_model = vae_adata.differential_expression(dge_data_sub, groupby = 'predictions')
    # get thresholds for markers
    dge_threshs = Toplevel()
    # lfc min
    lfc_min = StringVar()
    lfc_min.set("0.15")
    Message(dge_threshs, text = 'Minimum log-fold change:', width=250).pack()
    lfc_entry = Entry(dge_threshs, textvariable = lfc_min)
    lfc_entry.pack()
    # bayes min
    bayes_min = StringVar()
    bayes_min.set("2.5")
    Message(dge_threshs, text = 'Minimum bayes factor:', width=250).pack()
    bayes_entry = Entry(dge_threshs, textvariable = bayes_min)
    bayes_entry.pack()
    # non_zeros min
    non_zero = StringVar()
    non_zero.set("0.2")
    Message(dge_threshs, text = 'Fraction non-zero counts [0.0-1.0]:', width=250).pack()
    non_zero_entry = Entry(dge_threshs, textvariable = non_zero)
    non_zero_entry.pack()
    # max no of markers
    markers_max = StringVar()
    markers_max.set("5")
    Message(dge_threshs, text = 'Maximum number of markers per cluster:', width=250).pack()
    markers_entry = Entry(dge_threshs, textvariable = markers_max)
    markers_entry.pack()
    # wait for entry....
    wait_var = IntVar()
    def activate():
        wait_var.set(1)
    submit = Button(dge_threshs, text = "Apply thresholds", command = lambda:[activate(), dge_threshs.destroy()])
    submit.pack()
    submit.wait_variable(wait_var)
    # get thresholds
    lfc_min = float(lfc_min.get())
    bayes_min = float(bayes_min.get())
    markers_max = int(markers_max.get())
    non_zero = float(non_zero.get())
    # make marker list
    markers = {}
    cats = pd.Series(dge_data_sub.obs['predictions'].values, dtype="category").cat.categories
    for i, c in enumerate(cats):
        cid = "{} vs Rest".format(c)
        cell_type_df = dge_model.loc[dge_model.comparison == cid]
        cell_type_df = cell_type_df[cell_type_df.lfc_mean > lfc_min]
        cell_type_df = cell_type_df[cell_type_df["bayes_factor"] > bayes_min]
        cell_type_df = cell_type_df[cell_type_df["non_zeros_proportion1"] > non_zero]
        markers[c] = cell_type_df.index.tolist()[:markers_max]
    sc.pl.matrixplot(dge_data_sub, markers, cmap="inferno", groupby="predictions", standard_scale='var', dendrogram=True)
    # make spreadsheet
    make_stats_window(x=dge_model)
    reset_dge3()

def dge_betw_conditions():
    dge_data = adata_copy.copy()
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
    try:
        dge_model = vae_adata_red.differential_expression(dge_data_sub, groupby = 'condition')
    except:
        dge_model = vae_adata.differential_expression(dge_data_sub, groupby = 'condition')
    # get thresholds for markers
    dge_threshs = Toplevel()
    # lfc min
    lfc_min = StringVar()
    Message(dge_threshs, text = 'Minimum log-fold change:', width=250).pack()
    lfc_entry = Entry(dge_threshs, textvariable = lfc_min)
    lfc_entry.pack()
    # bayes min
    bayes_min = StringVar()
    Message(dge_threshs, text = 'Minimum bayes factor:', width=250).pack()
    bayes_entry = Entry(dge_threshs, textvariable = bayes_min)
    bayes_entry.pack()
    # non_zeros min
    non_zero = StringVar()
    Message(dge_threshs, text = 'Fraction non-zero counts [0.0-1.0]:', width=250).pack()
    non_zero_entry = Entry(dge_threshs, textvariable = non_zero)
    non_zero_entry.pack()
    # max no of markers
    markers_max = StringVar()
    Message(dge_threshs, text = 'Maximum number of markers per cluster:', width=250).pack()
    markers_entry = Entry(dge_threshs, textvariable = markers_max)
    markers_entry.pack()
    # wait for entry....
    wait_var = IntVar()
    def activate():
        wait_var.set(1)
    submit = Button(dge_threshs, text = "Apply thresholds", command = lambda:[activate(), dge_threshs.destroy()])
    submit.pack()
    submit.wait_variable(wait_var)
    # get thresholds
    lfc_min = float(lfc_min.get())
    bayes_min = float(bayes_min.get())
    markers_max = int(markers_max.get())
    non_zero = float(non_zero.get())
    # make marker list
    markers = {}
    cats = pd.Series(dge_data_sub.obs['condition'].values, dtype="category").cat.categories
    for i, c in enumerate(cats):
        cid = "{} vs Rest".format(c)
        cell_type_df = dge_model.loc[dge_model.comparison == cid]
        cell_type_df = cell_type_df[cell_type_df.lfc_mean > lfc_min]
        cell_type_df = cell_type_df[cell_type_df["bayes_factor"] > bayes_min]
        cell_type_df = cell_type_df[cell_type_df["non_zeros_proportion1"] > non_zero]
        markers[c] = cell_type_df.index.tolist()[:markers_max]
    sc.pl.matrixplot(dge_data_sub, markers, cmap="inferno", groupby="condition", standard_scale='var', dendrogram=True, layer = 'scvi_norm')
    make_stats_window(x=dge_model)
    reset_dge2()

# destroy subwindow
def close(x):
    x.destroy()

# User interface for plots
def interface_plotting():
    interface_plotting_window = Toplevel()
    # get logo
    global plot_image
    if getattr(sys, 'frozen', False):
        plot_image = PIL.Image.open(os.path.join(sys._MEIPASS, "files/Plot_options.png"))
    else:
        plot_image = PIL.Image.open("files/Plot_options.png")

    plot_image = plot_image.resize((1200,390))
    plot_image =  PIL.ImageTk.PhotoImage(plot_image)
    Label(interface_plotting_window, image = plot_image).pack(side=TOP)
    interface_plotting_window.title('Which plot do you want to create?')
    # checkbuttons
    Button(interface_plotting_window , text = 'Stacked barchart', width = 30, height = 1, command = lambda:[stacked_bar_plot()]).pack(side=LEFT)
    Button(interface_plotting_window , text = 'Boxplot by celltype', width = 54, height = 1, command = lambda:[boxplot_by_celltype()]).pack(side=LEFT)
    Button(interface_plotting_window, text = 'Boxplot by celltype and condition', width = 56, height = 1, command = lambda:[boxplot_by_condition()]).pack(side=LEFT)
    heatmap_options = ['Default lineage/subset heatmap', 'Custom heatmap']
    OptionMenu(interface_plotting_window, heatmap_choice, *heatmap_options, command = lambda _:[plot_heatmap()]).pack(side=LEFT)

# User interface for differential gene expression
def interface_dge():
    interface_dge_window = Toplevel()
    # image for dge
    try:
        global dge_image
        if getattr(sys, 'frozen', False):
            dge_image = PIL.Image.open(os.path.join(sys._MEIPASS, "files/dge_options.png"))
        else:
            dge_image = PIL.Image.open("files/dge_options.png")

        dge_image =  PIL.ImageTk.PhotoImage(dge_image)
        Label(interface_dge_window, image = dge_image).pack(side=TOP)
        # Button for "Between cells"
        OptionMenu(interface_dge_window, subset_choice5, *options5, command = lambda _:[dge_betw_celltypes()]).pack(side=LEFT)
        # Dropdown for "Between condition"
        OptionMenu(interface_dge_window, subset_choice4, *options4, command = lambda _:[dge_betw_conditions()]).pack(side=LEFT)
    except:
        print(colored('These subsets are only defined for the mapping workflow. Choose "All cells"'))

def keep_all():
    subset_choice2.set('Do not subset')
    subset_choice.set('Do not subset')

def interface_scVI_workflow():
    interface_scVI_window = Toplevel()
    # get logo
    global scVI_image
    # get logo
    if getattr(sys, 'frozen', False):
        scVI_image = PIL.Image.open(os.path.join(sys._MEIPASS, "files/acsf_logo_basic_workflow.png"))
    else:
        scVI_image = PIL.Image.open("files/acsf_logo_basic_workflow.png")
    scVI_image = scVI_image.resize((1300,200))
    scVI_image =  PIL.ImageTk.PhotoImage(scVI_image)
    wl = Label(interface_scVI_window, image = scVI_image)
    wl.pack()
    # make buttons
    Button(interface_scVI_window, text = 'I. a. Choose 10X data', command = lambda:[get_path()], width = 19, height = 1).pack(side=LEFT)
    Button(interface_scVI_window, text = 'II. a. Choose condition.tsv', command = lambda:[load_condition_csv()], width = 19, height = 1).pack(side=LEFT)
    Button(interface_scVI_window, text = 'II. b. Choose study.tsv', command = lambda:[load_study_csv()], width = 19, height = 1).pack(side=LEFT)
    Button(interface_scVI_window, text = 'II. c. Choose idents.tsv', command = lambda:[load_ident_csv()], width = 19, height = 1).pack(side=LEFT)
    Button(interface_scVI_window, text= 'III. Run basic workflow', command = lambda:[scVI_workflow()], width = 19, height = 1).pack(side=LEFT)
    Button(interface_scVI_window, text= 'IV. a. Plot UMAP', command = lambda :[keep_all(), plot_umap(), recall_full_data()], width = 19, height = 1).pack(side=LEFT)
    Button(interface_scVI_window, text = 'IV. b. Differential expression', command = interface_dge).pack(side=LEFT)
    Button(interface_scVI_window, text = 'VI. Rename leiden clusters', command = lambda: [get_new_names(), map_new_names(), keep_all()]).pack(side=LEFT)

def get_new_names():
    clusters = adata_copy.obs['predictions'].value_counts().index.tolist()
    global string_vars
    string_vars = {}
    for i in range(0, len(clusters)):
        locals()['clust' + str(i)] = StringVar()
        string_vars[i] = locals()['clust' + str(i)]
    rename_window = Toplevel()
    Label(rename_window, text = "Define new names for clusters.\nThey will be available as 'predictions' in UMAP plots.").pack()
    for i in range(0, len(clusters)):
        print(adata_copy.obs['predictions'].value_counts().index.tolist()[i])
        text_output = 'New name for clust' + str(i)
        Message(rename_window, text = 'Type new name for cluster ' + str(i) + ": " + '"' + str(adata_copy.obs['predictions'].value_counts().index.tolist()[i]) + '"', width=350).pack()
        locals()['rename_', str(i)] = Entry(rename_window, textvariable = string_vars[i])
        locals()['rename_', str(i)].pack()
    submit = Button(rename_window, text = "Apply new names", command = lambda:[submit_var.set(1), rename_window.destroy()])
    submit_var = IntVar()
    submit.pack()
    submit.wait_variable(submit_var)

def map_new_names():
    # map to predictions
    print(adata_copy)
    print(string_vars)
    print('start_loop')
    clusters = adata_copy.obs['predictions'].value_counts().index.tolist()
    for i in range(0, len(clusters)):
        print(adata_copy)
        print(colored('Changing name of Cluster...', 'yellow'), colored(i, 'white'))
        adata_copy.obs['predictions'] = np.where(adata_copy.obs['predictions'] == clusters[i], string_vars[i].get(), adata_copy.obs['predictions'])
    print(adata_copy.obs['predictions'].value_counts())

# make variables to store user-provided values
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
misc_choice = StringVar()
misc_choice.set('Miscallaneous functions')
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

def choose_misc_func(self):
    if misc_choice.get() == "Prepare data":
        prep_data()
    elif misc_choice.get() == "scanpy/scVI workflow":
        interface_scVI_workflow()
    elif misc_choice.get() == "Enrichment analysis":
        enrich_score()

# make button to start loading scvi data
butt1 = Button(root, text = 'I. a. Choose 10X data', command = lambda:[get_path()], width = 19, height = 1)
butt2 = Button(root, text = 'II. a. Choose condition.tsv', command = lambda:[load_condition_csv()], width = 19, height = 1)
butt3 = Button(root, text = 'II. b. Choose study.tsv', command = lambda:[load_study_csv()], width = 19, height = 1)
butt4 = Button(root, text = 'II. c. Choose idents.tsv', command = lambda:[load_ident_csv()], width = 19, height = 1)
butt5 = Button(root, text= 'III. Run mapping', command = lambda:[map_to_query()], width = 19, height = 1)
options3 = ['Do not subset', 'CD4 T cells', 'CD8 T cells', 'Myeloid cells', 'NK cells', 'B cells']
butt6 = OptionMenu(root, subset_choice2,  *options3, command = lambda _:[plot_umap(), reset_plot_umap(), recall_full_data()])
options = ['Do not subset', 'CD4 T cells', 'CD8 T cells', 'Myeloid cells', 'NK cells', 'B cells']
butt7 = OptionMenu(root, subset_choice, *options, command = lambda _:[subset_fun(), reset_quant(), interface_plotting(), recall_full_data()])
options2 = ['Annotations', 'UMAP coordinates', '.h5ad']
butt8 = OptionMenu(root, data_to_write, *options2, command = lambda _:[write_data(), reset_write()])
butt9 = Button(root, text = 'IV. c. Differential expression', command = interface_dge)
global options4
options4 = ['All cells', 'CD4 T cells', 'CD8 T cells', 'Myeloid cells', 'NK cells', 'B cells', ' B activated', 'B IL4R+', 'B atypical', 'Plasmacells', 'pDCs', 'Monocytes', 'BAM MRC1+', 'BAM EMP3+', 'MG CX3CR1+', 'MG CCL2+', 'MG TREM2hi', 'mDCs CD1c+', 'mDCs AXL+SIGLEC6+', 'mDCs CLEC9A+', 'NK bright', 'NK dim', 'TR-NK', 'ILC', 'MAIT', 'gdT Vd2+', 'gdT Vd2-', 'CD8 CM', 'CD8 EM HLA-DRA+', 'CD8 EM CD160+', 'CD8 TRM ITGA1+', 'CD8 TRM ITGA1-', 'CD8 CTL', 'Tfh', 'Th17', 'Th2/Th22', 'Th1', 'Tregs', 'CCR5high Th17.1', 'CD4 TEMRA']
global options5
options5 = ['All cells', 'CD4 T cells', 'CD8 T cells', 'Myeloid cells', 'NK cells', 'B cells']
options6 = ['Prepare data', 'scanpy/scVI workflow', 'Enrichment analysis']
butt10 = OptionMenu(root, misc_choice, *options6, command = choose_misc_func)

# make var to store genes to plot
gene_var = StringVar(root)
order_var = StringVar(root)
hue_var = StringVar(root)
stats_adata_var = StringVar(root)

# reset buttons
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
butt1.pack(side =LEFT)
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
