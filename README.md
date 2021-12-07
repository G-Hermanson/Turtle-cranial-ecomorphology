# Annotated R scripts and data associated with Hermanson et al. [Paper title]
This repository includes:

- 3D landmark coordinates of extant and fossil turtle skulls, for both the "full landmark dataset" and "partial landmark dataset" (See article). These can be found in the "Landmarks_full_dataset" and "Landmarks_partial_dataset" folders, respectively, within this repository

- Phylogenetic trees files in Newick format, used for phylogenetic regressions of skull (i) size and (ii) shape:
- Pereira_tree_pruned_full.tre: Molecular-based tree from Pereira et al. (2017, Methods in Ecol. & Evo.), pruned to our "full landmark dataset" sample
- Pereira_tree_pruned_partial.tre: Molecular-based tree from Pereira et al. (2017, Methods in Ecol. & Evo.), pruned to our "partial landmark dataset" sample
- Evers_tree_pruned: time-scaled composite tree including extant and extinct turtles, based on the consensus of Evers et al. (2019, PeerJ)
- Sterli_tree_pruned: time-scaled composite tree including extant and extinct turtles, based on the MkA tree of Sterli et al. (2018, J. Vert. Pal.)
  
- R scripts used to run the analyses
- Script 1 (load 3D, GMM, Hypothesis tests - Full landmark dataset).R: used to load 3D landmark data of the "full landmark dataset", run geometric morphometric analyses (GPA, PCA), and run hypothesis tests (D-PGLS, pGLS, 2B-PLS)
- Script 2 (load 3D, GMM, Hypothesis tests, pFDA - Partial landmark dataset).R: used to load 3D landmark data of the "partial landmark dataset", run geometric morphometric analyses (GPA, PCA), run hypothesis tests (D-PGLS), and predict traits for fossil turtles using phylogenetic flexible discriminant analyses (pFDA; Motani & Schmitz 2011, Evolution)
- Script 3 (plots for graphical visualisation of results).R: used for graphical visualisation of results, including (i) ordination plots (PCA, D-PGLS, pGLS, 2B-PLS), (ii) shape deformation 3D plots (PCA, D-PGLS, 2B-PLS), and (iii) trait-mapping on trees
- Custom R functions.R: this script contains (i) custom R code used to retrieve results from D-PGLS objects, and (ii) custom functions from Motani and Schmitz (2011, Evolution) used to perform pFDA
