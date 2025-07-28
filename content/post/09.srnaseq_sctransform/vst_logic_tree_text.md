# 🔧 `vst()` Function Logic Tree with Triggers (Seurat/scTransform-style)

## 1. Preprocessing and Argument Setup

**Preprocessing Flow (Text Tree):**
```
Start vst() function
├── Check vst.flavor
│   ├── If vst.flavor == "v2":
│   │   ├── method ← "glmGamPoi_offset"
│   │   ├── exclude_poisson ← TRUE
│   │   ├── If min_variance == -Inf:
│   │   │   └── min_variance ← "umi_median"
│   │   └── If n_cells is NULL:
│   │       └── n_cells ← 2000
│   │
│   └── If vst.flavor == "seurat_v3":
│       ├── method ← "nb_fast"
│       ├── exclude_poisson ← FALSE
│       └── Leave min_variance, n_cells unchanged
│
├── Handle deprecated arguments:
│   ├── If verbose is provided:
│   │   └── show_progress ← verbose
│   └── Keep show_progress for downstream
│
└── Dependency check:
    └── If method in ("glmGamPoi", "glmGamPoi_offset"):
        └── Check for installed glmGamPoi package
```

**Decision Tree Logic:**
```yaml
preprocessing:
  vst_flavor_v2:
    method: "glmGamPoi_offset"
    exclude_poisson: true
    min_variance_default: "umi_median"  # if -Inf
    n_cells_default: 2000  # if NULL
    
  vst_flavor_seurat_v3:
    method: "nb_fast"
    exclude_poisson: false
    # Keep existing min_variance, n_cells
    
  deprecated_handling:
    verbose: show_progress  # Map old arg to new
    
  dependency_check:
    glmGamPoi_methods: ["glmGamPoi", "glmGamPoi_offset"]
    required_package: "glmGamPoi"
```

## 2. Cell and Gene Filtering

**Cell and Gene Filtering Flow (Text Tree):**
```
Start Filtering Process
├── Cell Filtering:
│   ├── Remove cells with total UMI < min_cells_total
│   └── If n_cells is not NULL:
│       └── Randomly subsample n_cells cells
│
├── Gene Filtering (Sequential):
│   ├── If min_cells is set:
│   │   └── Keep genes expressed in ≥ min_cells cells
│   ├── If min_expr is set:
│   │   └── Keep genes with total counts ≥ min_expr
│   └── If max_cells is set:
│       └── Remove genes expressed in > max_cells cells
│
└── Return filtered UMI matrix
```

**Filtering Parameters:**
```r
# Cell filtering
filter_cells <- function(umi_matrix, min_cells_total, n_cells = NULL) {
  # Remove low-UMI cells
  cell_totals <- colSums(umi_matrix)
  keep_cells <- cell_totals >= min_cells_total
  
  # Optional subsampling
  if (!is.null(n_cells) && sum(keep_cells) > n_cells) {
    keep_indices <- sample(which(keep_cells), n_cells)
    keep_cells <- rep(FALSE, ncol(umi_matrix))
    keep_cells[keep_indices] <- TRUE
  }
  
  return(umi_matrix[, keep_cells])
}

# Gene filtering
filter_genes <- function(umi_matrix, min_cells = NULL, min_expr = NULL, max_cells = NULL) {
  keep_genes <- rep(TRUE, nrow(umi_matrix))
  
  if (!is.null(min_cells)) {
    genes_per_cell <- rowSums(umi_matrix > 0)
    keep_genes <- keep_genes & (genes_per_cell >= min_cells)
  }
  
  if (!is.null(min_expr)) {
    gene_totals <- rowSums(umi_matrix)
    keep_genes <- keep_genes & (gene_totals >= min_expr)
  }
  
  if (!is.null(max_cells)) {
    genes_per_cell <- rowSums(umi_matrix > 0)
    keep_genes <- keep_genes & (genes_per_cell <= max_cells)
  }
  
  return(umi_matrix[keep_genes, ])
}
```

## 3. Model Initialization (Step 1 Genes)

**Model Initialization Flow (Text Tree):**
```
Choose genes_step1
├── If n_genes is set:
│   └── Randomly choose n_genes from filtered genes
└── Else:
    └── Use all filtered genes

Model Parameter Setup
├── If fix_slope OR fix_intercept:
│   ├── Compute gene_mean ← rowMeans(UMI)
│   ├── Compute mean_cell_sum ← mean(colSums(UMI))
│   ├── If use_geometric_mean_offset:
│   │   └── Replace gene_mean with row_gmean(UMI)
│   └── Create model_pars_fixed:
│       ├── intercept ← NA (estimated by offset)
│       ├── slope ← log(gene_mean) - log(mean_cell_sum)
│       └── log_overdispersion ← log(10)
│
└── Fit initial model for genes_step1:
    ├── If method == "glmGamPoi_offset":
    │   └── Fit overdispersion via glmGamPoi::overdispersion_mle()
    └── Else:
        └── Use default NB model (glm.nb or internal fitting)
```

**Model Parameter Setup:**
```r
initialize_model <- function(umi_matrix, fix_slope = FALSE, fix_intercept = FALSE,
                           use_geometric_mean_offset = FALSE, n_genes = NULL) {
  
  # Gene selection for step 1
  if (!is.null(n_genes) && n_genes < nrow(umi_matrix)) {
    genes_step1 <- sample(nrow(umi_matrix), n_genes)
  } else {
    genes_step1 <- seq_len(nrow(umi_matrix))
  }
  
  # Fixed parameter setup
  if (fix_slope || fix_intercept) {
    if (use_geometric_mean_offset) {
      gene_mean <- apply(umi_matrix, 1, function(x) exp(mean(log(x[x > 0]))))
    } else {
      gene_mean <- rowMeans(umi_matrix)
    }
    
    mean_cell_sum <- mean(colSums(umi_matrix))
    
    # [intercept, slope, log_overdispersion]
    # Intercept estimated by offset; slope fixed at log(10)
    model_pars_fixed <- cbind(
      intercept = NA,
      slope = log(gene_mean) - log(mean_cell_sum),
      log_overdispersion = log(10)
    )
    
    return(list(
      genes_step1 = genes_step1,
      model_pars_fixed = model_pars_fixed
    ))
  }
  
  return(list(genes_step1 = genes_step1))
}
```

## 4. Model Fit for Remaining Genes (Binned)

**Binned Gene Fitting Flow (Text Tree):**
```
Bin Remaining Genes
├── Split genes into batches (~1000 genes per bin)
└── For each bin:
    ├── Reuse regressor matrix from Step 1
    ├── Fitting strategy:
    │   ├── If fix_slope OR fix_intercept:
    │   │   ├── Use model_pars_fixed from Step 1
    │   │   └── Predict: mu_bin ← exp(offset + slope × log_umi)
    │   │
    │   └── Else (fit from scratch):
    │       ├── If method == "glmGamPoi_offset":
    │       │   └── Fit via glmGamPoi with offset
    │       ├── If method == "nb_fast":
    │       │   └── Fit NB model via glm.nb()
    │       └── Else:
    │           └── Use custom NB fitting routine
    │
    └── Store bin results (θ, intercept, slope per gene)
```

**Binned Fitting Implementation:**
```r
fit_genes_binned <- function(umi_matrix, regressor_matrix, model_pars_fixed = NULL,
                           method = "glmGamPoi_offset", bin_size = 1000) {
  
  n_genes <- nrow(umi_matrix)
  gene_bins <- split(seq_len(n_genes), ceiling(seq_len(n_genes) / bin_size))
  
  results <- vector("list", length(gene_bins))
  
  for (i in seq_along(gene_bins)) {
    bin_genes <- gene_bins[[i]]
    bin_data <- umi_matrix[bin_genes, , drop = FALSE]
    
    if (!is.null(model_pars_fixed)) {
      # Use fixed parameters
      results[[i]] <- fit_with_fixed_params(bin_data, regressor_matrix, 
                                          model_pars_fixed[bin_genes, ])
    } else {
      # Fit from scratch
      results[[i]] <- switch(method,
        "glmGamPoi_offset" = fit_glmGamPoi_offset(bin_data, regressor_matrix),
        "nb_fast" = fit_nb_fast(bin_data, regressor_matrix),
        fit_default_nb(bin_data, regressor_matrix)
      )
    }
  }
  
  return(do.call(rbind, results))
}
```

## 5. Variance Stabilization and Residual Calculation

**Variance Stabilization and Residual Calculation Flow (Text Tree):**
```
For each gene and cell pair (g,c):
├── Compute expected value: μ̂_gc from fitted model
├── Compute variance estimate: V̂_gc = μ̂ + μ̂²/θ
├── Handle min_variance:
│   ├── If min_variance == 'umi_median':
│   │   └── min_var ← (median(nonzero_UMIs) / 5)²
│   ├── If numeric value provided:
│   │   └── min_var ← provided value  
│   └── Else:
│       └── min_var ← 0
│
├── Calculate Pearson residual:
│   └── residual_gc ← (y_gc - μ̂_gc) / sqrt(max(V̂_gc, min_var))
│
└── Cap extreme residuals:
    ├── If residual > 50: set to 50
    └── If residual < -50: set to -50
```

**Residual Calculation:**
```r
calculate_pearson_residuals <- function(umi_matrix, mu_matrix, theta_vector, 
                                      min_variance = NULL) {
  
  # Handle min_variance settings
  if (is.character(min_variance) && min_variance == "umi_median") {
    nonzero_umis <- umi_matrix[umi_matrix > 0]
    min_var <- (median(nonzero_umis) / 5)^2
  } else if (is.numeric(min_variance)) {
    min_var <- min_variance
  } else {
    min_var <- 0
  }
  
  # Calculate variance for each gene-cell pair
  variance_matrix <- mu_matrix + (mu_matrix^2) / theta_vector
  variance_matrix <- pmax(variance_matrix, min_var)
  
  # Pearson residuals
  residuals <- (umi_matrix - mu_matrix) / sqrt(variance_matrix)
  
  # Cap extreme residuals
  residuals[residuals > 50] <- 50
  residuals[residuals < -50] <- -50
  
  return(residuals)
}
```

## 6. Output Formatting

**Output Formatting Flow (Text Tree):**
```
Collect and Format Results
├── Main outputs:
│   ├── Pearson residual matrix (genes × cells)
│   ├── Model parameters per gene:
│   │   ├── θ (overdispersion parameter)
│   │   ├── intercept 
│   │   └── slope
│   └── Fitted values matrix: μ (genes × cells)
│
├── Metadata and diagnostics:
│   ├── Gene/cell filtering information:
│   │   ├── Number of genes removed
│   │   ├── Number of cells removed  
│   │   ├── Final gene count
│   │   └── Final cell count
│   │
│   └── Optional diagnostics:
│       ├── Overdispersion metrics
│       ├── Model convergence info
│       └── Runtime statistics
│
└── Return structured vst_result object
```

**Output Structure:**
```r
format_vst_output <- function(residuals, model_params, mu_matrix, 
                             filter_info, diagnostics = NULL) {
  
  vst_output <- list(
    # Main results
    residuals = residuals,  # genes x cells matrix
    
    # Model parameters
    model_pars = data.frame(
      gene = rownames(residuals),
      theta = model_params$theta,
      intercept = model_params$intercept,
      slope = model_params$slope
    ),
    
    # Fitted values
    mu = mu_matrix,
    
    # Filtering information
    filter_info = list(
      genes_removed = filter_info$genes_removed,
      cells_removed = filter_info$cells_removed,
      n_genes_final = nrow(residuals),
      n_cells_final = ncol(residuals)
    ),
    
    # Optional diagnostics
    diagnostics = diagnostics
  )
  
  class(vst_output) <- "vst_result"
  return(vst_output)
}
```

## ✅ Summary of Key Triggers

| **Condition** | **Result** |
|---------------|------------|
| `vst.flavor == "v2"` | Switch to `glmGamPoi_offset` + modern defaults |
| `fix_intercept` or `fix_slope` | Construct offset-based model instead of full fitting |
| `method == "glmGamPoi_offset"` | Call `glmGamPoi::overdispersion_mle()` for θ |
| `min_variance == 'umi_median'` | Set variance floor so Pearson residual ≤ 5 for 1 UMI genes |

```yaml
trigger_summary:
  flavor_triggers:
    v2: 
      method: "glmGamPoi_offset"
      exclude_poisson: true
      defaults: ["umi_median", 2000]
    seurat_v3:
      method: "nb_fast" 
      exclude_poisson: false
      
  modeling_triggers:
    fixed_parameters: "offset_based_fitting"
    glmGamPoi_offset: "overdispersion_mle"
    umi_median: "variance_floor_protection"
    
  output_triggers:
    residual_capping: [-50, 50]
    diagnostics: "optional_metrics"
```