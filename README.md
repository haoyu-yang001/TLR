# TLRmediation

`TLRmediation` is an R package for large-scale causal mediation testing under a composite null hypothesis, motivated by epigenome-wide studies with a large number of candidate mediators.

The package implements the **tail likelihood ratio (TLR)** procedure for testing whether an exposure affects an outcome through a mediator, together with supporting utilities for:

- composite-null proportion estimation
- tail-parameter estimation
- empirical-null correction
- comparison with existing methods such as **Sobel**, **DACT**, and **HDMT**
- simulation, visualization, and performance evaluation

---

## Background

In large-scale mediation analysis, each candidate mediator is associated with two marginal tests:

- an **exposure--mediator** association test
- a **mediator--outcome** association test

The null hypothesis of no mediation effect is therefore **composite**, since mediation is absent whenever at least one of the two component effects is zero. This creates major challenges for valid and powerful large-scale inference.

The main goal of `TLRmediation` is to provide a practical implementation of the **tail likelihood ratio (TLR)** framework for large-scale causal mediation testing, while also including several existing methods for comparison.

---

## Main functionality

### Core TLR procedures

- `null_estimation()`  
  Estimate composite-null mixture proportions from paired p-values.

- `alphas_estimate_tail()`  
  Estimate tail approximation parameters for alternative p-value distributions.

- `TLR()`  
  Apply the TLR multiple-testing procedure and return selected hypotheses.

- `TLR_thre()`  
  Thresholding step for the TLR procedure.

- `TLR_pvalues()` / `TLR_pvalues_trimmed()`  
  Compute TLR-based p-values under estimated null distributions.

### Comparison methods

- `Sobel()`  
  Classical Sobel mediation test.

- `DACT_liu()`  
  DACT-style combined p-values.

- `MDACT()` / `MDACT_pvalues()`  
  MDACT-style procedures.

- `HDMT_emp()` / `HDMT_asy()` / `HDMT_control()`  
  Empirical and asymptotic HDMT procedures.

### Empirical-null and auxiliary estimation

- `JCCorrect()`  
  Jin--Cai empirical-null correction.

- `EfronCorrect()`  
  Efron empirical-null correction.

- `nullParaEst()` / `nonnullPropEst()`  
  Empirical-null location/scale and non-null proportion estimation.

- `JC_prop_estmate()`  
  Estimate joint mixture proportions from paired Z-statistics.

### Simulation and evaluation

- `generate_pvalues_indp_from_norm()`  
  Simulate paired Z-statistics and p-values under a four-component mixture model.

- `index_fpn_fdr_power_func()`  
  Summarize FDR and power.

- `size_func()`  
  Summarize size and power.

### Plotting

- `qqunif.plot()`  
  QQ plot against the uniform null distribution.

---

## Installation

During development, you can load the package locally with:

```r
devtools::load_all()
