from stlearn import Read10X, ReadXenium, ReadMERFISH, ReadSeqFish, ReadSlideSeq
import os
import scanpy as sc
import pandas as pd
import json as js

def read_10x(self, data_type, directory, xenium_imgpath = None):
    exception_occurred = False
    try:
        if os.path.isdir(directory):
            if data_type == '10X Visium':
                self.parent.adata = Read10X(directory)
            elif data_type == '10X Xenium':
                self.parent.adata = ReadXenium(feature_cell_matrix_file= directory + "\\cell_feature_matrix.h5",
                cell_summary_file= directory + "\\cells.csv.gz",
                image_path=xenium_imgpath,
                library_id= os.path.basename(directory),
                scale=1,
                spot_diameter_fullres=15
                )
                # cells_df = df = pd.read_csv(directory + "\\cells.csv.gz")
                # cells_df.set_index(self.parent.adata.obs_names, inplace=True)
                # self.parent.adata.obs = df.copy()
                # self.parent.adata.obsm["spatial"] = self.parent.adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()
            elif data_type == 'MERFISH':
                self.parent.adata = ReadMERFISH(count_matrix_file=directory+"\\count_matrix.csv", spatial_file=directory+"\\spatial_file.xlsx")
            elif data_type == 'SeqFISH':
                self.parent.adata = ReadSeqFish(count_matrix_file=directory+"\\count_matrix.matrix", spatial_file=directory+"\\spatial_file.csv")
                self.parent.adata.var.index.name = None
            elif data_type == 'SlideSeq':
                self.parent.adata = ReadSlideSeq(count_matrix_file=directory+"\\count_matrix.count", spatial_file=directory+"\\spatial_file.idx")
            if 'in_tissue' not in self.parent.adata.obs.columns:
                self.parent.adata.obs['in_tissue'] = 1
            if 'raw_counts' not in self.parent.adata.layers:
                self.parent.adata.layers['raw_counts'] = self.parent.adata.X.copy()
                print("Copied counts to layer: [raw_counts]")

        else:
            self.parent.adata = sc.read_h5ad(directory)
            if 'spatial' not in self.parent.adata.uns or 'spatial' not in self.parent.adata.obsm:
                self.parent.adata = None
                raise Exception("The key 'spatial' is missing in either uns or obsm.")
            
        setInitialVars(self)

    except Exception as e:
        print(e)
        exception_occurred = True
        self.root.command_queue.put(('destroy_progress_window', {'value':'ind', 'title':'Failed to import!', 'msg':'The dataset has failed to import. \n For details on possible causes, please consult the documentation.', 'show':'fail'}))
    finally:
        if not exception_occurred:
            self.root.command_queue.put(('destroy_progress_window', {'value':'ind', 'title':'Successfully imported!', 'msg':'The dataset has been successfully imported.', 'show':'success'}))
            self.root.command_queue.put(('managedata_refreshstats', ()))

def setInitialVars(self):
    self.parent.gene_names = self.parent.adata.var.index.tolist()
    self.parent.default_total_spots.set(len(self.parent.adata.obs))
    self.parent.total_spots.set(len(self.parent.adata.obs))
    self.parent.default_total_genes.set(len(self.parent.adata.var))
    self.parent.total_genes.set(len(self.parent.adata.var))
    self.parent.dataname = list(self.parent.adata.uns['spatial'].keys())[0]                                             #Name of the dataset
    self.parent.imgsize = self.parent.adata.uns['spatial'][self.parent.dataname]['images']['hires'].shape[:2][::-1]     #Set the size of the original hires image in the dataset

    
    # self.parent.adata.var_names_make_unique()
    # self.parent.adata.layers["counts"] = self.parent.adata.X.copy().tocsr()                                             #Save the raw counts
    human_mt = self.parent.adata.var_names.str.startswith("MT-")
    mouse_mt = self.parent.adata.var_names.str.startswith("mt-")
    self.parent.adata.var["mt"] =  human_mt if human_mt.sum() > mouse_mt.sum() else mouse_mt
    sc.pp.calculate_qc_metrics(self.parent.adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    self.parent.ngenesbycounts.set(self.parent.adata.obs["n_genes_by_counts"].max())
    self.parent.total_counts.set(self.parent.adata.obs["total_counts"].max())
    self.parent.mit_counts.set(self.parent.adata.obs["pct_counts_mt"].max())
    self.parent.dataset_displayname.set(self.parent.dataname)

def saveas(self, path):
    exception_occurred = False
    try:
        self.parent.adata.write_h5ad(path + '\\adata.h5ad')
    except Exception as e:
        print('Exception in saving anndata: ', e)
        exception_occurred = True
        self.root.command_queue.put(('destroy_progress_window', {'value':'ind', 'title':'Failed to Export', 'msg':'The Anndata has failed to export.', 'show':'fail'}))
    finally:
        if not exception_occurred:
            self.root.command_queue.put(('destroy_progress_window', {'value':'ind','title':'Successful', 'msg':'Successfully exported the selected files.', 'show':'success'}))