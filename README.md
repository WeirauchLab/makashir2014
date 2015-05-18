README meta_analysis_code.R (makashir2014)
==========================================

The `meta_analysis_code.R` program takes m gene expression datasets from _m_
different studies as an input and performs meta-analysis of differential gene
co-expression on pairwise combinations of genes.

All the input data should be located in one single folder; the path of this
folder should be specified as `data.loc` in the code. Similarly, the output
path for results should be specified as `res.loc` in the code.

Datasets must be in tab separated matrix format, each row containing
expression of that gene and each column indicating that sample. 

There must be the only two unique column names in a dataset (each referring to
one condition e.g. cancer, healthy) and these should be the same across all
datasets. The first column must contain gene names.

There must be at least 4 samples in each study.

The following output files will be produced: 

1. `gene list`

   contains a list of all genes used for meta-analysis

2. `gene-pairs`

   contains the list of all gene pairs on which differential co-expression
   analysis is performed

3. `counts`

   matrix indicating the number of studies / datasets involved in
   differential co-expression meta-analysis of that gene pairs

4. `q_scores`

   matrix indicating the q scores (meta-analysis standard normal scores) for
   that gene pair

5. `p_values`
     matrix indicating the p values of  differential co-expression
     meta-analysis for that gene pair 

Note that this program does not adjust the p-values for multiple hypothesis
testing, this needs to be done separately. 
