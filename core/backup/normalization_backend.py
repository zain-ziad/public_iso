import scanpy as sc
import stlearn as st
from core.assets.themecolors import COLORS
from stlearn.image_preprocessing.model_zoo import encode, Model
from typing import Optional, Union
from anndata import AnnData
import numpy as np
from stlearn._compat import Literal
from PIL import Image
import pandas as pd
from customtkinter import IntVar
from tqdm import tqdm
from math import ceil

_CNN_BASE = Literal["resnet50", "vgg16", "inception_v3", "xception"]


def NormalizeTotal(self, adata, action_left):
    try:
        sc.pp.normalize_total(adata, inplace=True)
        self.parent.isNormalized.set(True)
    finally:
        if not action_left:
            self.root.command_queue.put(('destroy_progress_window', ('Successful!', 'Normalization done.', 'success')))  

def Log1p_Transform(self, adata, action_left):
    try:
        sc.pp.log1p(adata)
        self.parent.isLog1p.set(True)
    finally:
        self.log1p_button.action_completed()
        if not action_left:
            self.root.command_queue.put(('destroy_progress_window', ('Successful!', 'Log1p transformation done.', 'success'))) 

def Scale_Transform(self, adata, action_left):
    try:
        sc.pp.scale(adata)
        self.parent.isScaled.set(True)
    finally:
        self.scale_button.action_completed()
        if not action_left:
            self.root.command_queue.put(('destroy_progress_window', ('Successful!', 'Data successfully scaled.', 'success')))

def run_stSMEnormalization(self, adata, action_left):

    def ExtractFeature_Progress(
        adata: AnnData,
        cnn_base: _CNN_BASE = "resnet50",
        n_components: int = 50,
        verbose: bool = False,
        copy: bool = False,
        seeds: int = 1
    ) -> Optional[AnnData]:
        
        feature_dfs = []
        model = Model(cnn_base)

        if "tile_path" not in adata.obs:
            raise ValueError("Please run the function stlearn.pp.tiling")

        #   To properly increment the progress, we will use the following variables #
        end = len(adata)
        current_value = 0
        step_size = end / 100.0

        with tqdm(
            total=len(adata),
            desc="Extract feature",
            bar_format="{l_bar}{bar} [ time left: {remaining} ]",
        ) as pbar:
            for i, (spot, tile_path) in enumerate(adata.obs["tile_path"].items()):
                tile = Image.open(tile_path)
                tile = np.asarray(tile, dtype="int32")
                tile = tile.astype(np.float32)
                tile = np.stack([tile])
                if verbose:
                    print("extract feature for spot: {}".format(str(spot)))
                features = encode(tile, model)
                feature_dfs.append(pd.DataFrame(features, columns=[spot]))
                pbar.update(1)
                new_value = int(i / step_size)
                if new_value > current_value:
                    # Calculate time remaining and percentage
                    elapsed = pbar.format_dict['elapsed']
                    rate = pbar.format_dict['rate']
                    remaining = (pbar.total - pbar.n) / rate if rate else float('inf')
                    percentage = int((pbar.n / pbar.total) * 100)

                    #Update progress
                    current_value = new_value
                    self.root.detprogressmodal.update_det_progress(mins=ceil(remaining / 60), percent=percentage, elapsed = int(elapsed//60))

        feature_df = pd.concat(feature_dfs, axis=1)
        adata.obsm["X_tile_feature"] = feature_df.transpose().to_numpy()
        from sklearn.decomposition import PCA
        pca = PCA(n_components=n_components, random_state=seeds)
        pca.fit(feature_df.transpose().to_numpy())
        adata.obsm["X_morphology"] = pca.transform(feature_df.transpose().to_numpy())

        return adata if copy else None
    
    try:
        st.pp.tiling(adata)
        ExtractFeature_Progress(adata)
        st.em.run_pca(adata,n_comps=50)
        st.spatial.SME.SME_normalize(adata, use_data="raw")
        adata.X = adata.obsm['raw_SME_normalized']
        self.parent.isNormalized.set(True)
    finally:
        if not action_left:
            self.root.command_queue.put(('destroy_det_progress_window', ('Successful!', 'Successfully completed!', 'success')))
        else:
            self.root.command_queue.put(('destroy_det_progress_window', {'releaseAppWindow':False}))

def run_stSMEnormalization_noFeatureExtract(self, adata, action_left):
    try:
        st.spatial.SME.SME_normalize(adata, use_data="raw")                                                                                                      #Else, just apply stSME
        adata.X = adata.obsm['raw_SME_normalized']
        self.parent.isNormalized.set(True)
    finally:
        if not action_left:
            self.root.command_queue.put(('destroy_det_progress_window', ('Successful!', 'stSME normalization done.', 'success')))