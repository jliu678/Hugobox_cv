---
title: 🧬 Dynamic RNA velocity model-- (8) scVelo pipeline  
# draft: True
summary: Here delves into. 
date: 2025-05-28
authors:
  - admin
tags:
  - scRNAseq RNA Velocity, 
  - Dynamic model
  - scVelo
  - 
image:
  caption: 'Image credit: [**Logan Voss on Unsplash**](https://unsplash.com)'
---

The dynamic RNA velocity model makes assumptions that are not always satisfied in real-world data, which warns that default pipelines may fail to capture true RNA velocity. Our previous six blog posts on effectively applying this model have equipped us to identify and thoughtfully fine-tune key (hyper)parameters and processing steps in scVelo, allowing for more accurate and meaningful analyses.

However, to rigorously validate these adjustments, we need high-quality simulated data with known ground truth. That’s why this seventh blog delves into the Gillespie Stochastic Simulation Algorithm (SSA)-- a Monte Carlo method well suitable to simulate biochemical reaction systems where molecule numbers are low and stochastic effects are significant.

## Run scVelo on simulated data

### Wrangling the simulated data

Simulated data were downloaded [here](https://github.com/Spencerfar/LatentVelo/blob/main/synthetic_datasets/bifurcation_2batches.h5ad)

```python
# Import necessary libraries
import scvelo as scv
import scanpy as sc
import logging
import sys, re, random
import numpy as np
from datetime import datetime
from scipy import sparse

b2b=sc.read('/data/latentVelo_github/bifurcation_2batches.h5ad')
# Below verify the b2b.X is spliced
np.allclose((b2b.X - b2b.layers['counts_spliced']).data, 0)
```

I decide to:
- store raw data in a layer
- store spliced counts in both adata.X and adata.layers
- store unspliced counts in adata.layers

It is according to the structure of scVelo object, said [here](https://github.com/theislab/scvelo/issues/525#issuecomment-891056397)
> The current workflow works with the spliced counts in adata.X and expects the layers 'unspliced' and 'spliced'. You can, of course, store the raw count matrix in adata.layers as well if you need it.

and [scvelo stores unspliced and spliced counts (spliced both in adata.X and adata.layers)](https://github.com/scverse/scanpy/issues/1860#issuecomment-873990253)

```python
b2b.layers['spliced'] = b2b.layers['counts_spliced']
b2b.layers['unspliced'] =b2b.layers['counts_unspliced']
```

### Verify `milestone` and its relation to `traj_progression`
```python
b2b.uns['traj_progressions']['edge'] = (
    b2b.uns['traj_progressions']['from'] + '->' + b2b.uns['traj_progressions']['to']
)

b2b.uns['traj_progressions']['edge'].index = [
    f"cell{int(i) + 1}" for i in b2b.uns['traj_progressions']['edge'].index
]

b2b.obs['edge']=b2b.uns['traj_progressions']['edge'].loc[b2b.obs['milestone'].index]

import pandas as pd
pd.crosstab(b2b.obs['milestone'], b2b.obs['edge']).T.to_csv(
    'milestone_edge_crosstab.csv', 
                       index=True,          
                       header=True,          
                       index_label='edge') 
```
The result below shows the ground truth of the milestone trajectories, which contain two branches-- A->B->D and A->C->E.

{{< table path="milestone_edge_crosstab.csv" header="true" caption="" class="table-striped" >}}

### Default scVelo pipeline fail to capture the ground truth
```python
logging.info('1c.Preprocessing,PCA,findNeighbour,Calculate Moments')
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

logging.info('1d.EM fitting and calculate scVelo dynamic RNA velocity')
scv.tl.recover_dynamics(adata,n_jobs=81)
scv.tl.velocity(adata, mode='dynamical')

logging.info('1e.calculate Markov chain of the differentiation process')
scv.tl.velocity_graph(adata)

logging.info("2.velocity embeding(independent of root cells)")
scv.pl.velocity_embedding_stream(adata, basis='umap',color='milestone',
                                 legend_loc='right margin',show=True,
                                 save='b2b_default_embedding.png')

logging.info('2a.infer root cell, and calculate gene-shared latent time corrected by neighborhood convolution')
scv.tl.latent_time(adata)# identical to recover_latent_time introduced in scVelo paper
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80,dpi=300, show=True,
               save='b2b_default_latent-time.png')
```

{{< figure
  src="scvelo_b2b_ground-truth_default-scvelo.png"
  alt=""
  link=""
  caption="Default scVelo pipeline fail to capture the ground truth"
  class="ma0 w-75"
>}}

### scVelo performance depends on dataset 'complexity' and embedding space
When embeded in PCA space:
{{< figure
  src="b2b_velocity-stream_pca.png"
  alt=""
  link=""
  caption="scVelo performance on datasets with varied 'complexity' in PCA space"
  class="ma0 w-75"
>}}

When embeded in UMAP space:
{{< figure
  src="b2b_velocity-stream_umap.png"
  alt=""
  link=""
  caption="scVelo performance on datasets with varied 'complexity' in UMAP space"
  class="ma0 w-75"
>}}

### Prior knowledge of root cells improves accuracy of latent time
{{< figure
  src="b2b_latent-time.png"
  alt=""
  link=""
  caption="scVelo performance on datasets with varied 'complexity' in UMAP space"
  class="ma0 w-75"
>}}


### Borrow velocity, PCA and/or neighbor from larger population of cells

{{< figure
  src="b2b_borrow3.png"
  alt=""
  link=""
  caption=""
  class="ma0"
>}}

## Observations and discussions

1. The default scVelo pipeline fails to capture the ground truth of the bifurcation process, as shown in the figures above.

Prior knowledge of root cells and ... is required. Otherwise, no way to evaluate the accuracy of the inferred velocity stream or latent time.

2. The performance of scVelo is highly dependent on the dataset's 'complexity' and the embedding space used. The results show that scVelo performs better in PCA space than in UMAP space, which may be due to the fact that PCA distorted data less than UMAP. But we have to consider that UMAP often resolve high dimensional data better than PCA as PCA is linear, UMAP is nonlinear.
- High-dimensional biological data (e.g., single-cell RNA-seq) often lies on nonlinear manifolds. UMAP can model this nonlinear structure, such as branching trajectories, subtle gradients and clusters that PCA would flatten or blur.

3. Borrowing velocity, PCA, and/or neighbor information from a larger population of cells can significantly improve the performance of scVelo, as shown in the last figure. This suggests that scVelo can benefit from leveraging information from other populations to enhance its analysis.

## real-world data: tumor development
{{< figure
  src="ut.png"
  alt=""
  link=""
  caption=""
  class="ma0"
>}}
