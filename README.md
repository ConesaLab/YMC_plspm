# YMC_plspm

This is the code that was used to preprocess and analyze the data corresponding to the paper "Metabolite accumulation mediates the shift between the High Oxygen Consumption
and Low Oxygen Consumption phases in the Yeast Metabolic Cycle".

In this article we statistically integrated Transcripts, Chromatin and Metabolites in a single model in the Yeast Metabolic Cycle. Successfully modeling this dynamic system, and defining subsets of metabolites that impact the cycle progression.

The different scripts were used for the following purposes:

Splines_matrices.R -> Generate the spline matrices from the processed data for each of the omic datasets. Matching the timepoints.
Clustering_plots.R -> Cluster each omic separately and weight the optimal number of clusters (Figure 2)
PLS_splines.R -> Partial Least Squares Modeling of Metabolites and Gene Expression (Figure 2)
prepare_plspm_datasets.R -> Create the PLS-PM datasets and prepare the model (Figure 4)
plspm_models_221120.R -> Run the PLSPM model (Figure 4)
monte_carlo_plspm.R -> Run the iterative plspm process (Figure 4)
early_late_metabolites.R -> Plot the metabolite heatmaps (Figure 5)
HM_correlations.R -> Address if the H3K18ac and H3K9ac signals are correlated
crossloadings_assessment.R -> Crossloading of Figure S4
