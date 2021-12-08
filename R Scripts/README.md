# R scripts
This folder contains R scripts used to run the analyses

- Script 1 (load 3D, GMM, Hypothesis tests - Full landmark dataset).R: used to load 3D landmark data of the "full landmark dataset", run geometric morphometric analyses (GPA, PCA), and run hypothesis tests (D-PGLS, pGLS, 2B-PLS);
- Script 2 (load 3D, GMM, Hypothesis tests, pFDA - Partial landmark dataset).R: used to load 3D landmark data of the "partial landmark dataset", run geometric morphometric analyses (GPA, PCA), run hypothesis tests (D-PGLS), and predict traits for fossil turtles using phylogenetic flexible discriminant analyses (pFDA; Motani & Schmitz 2011, Evolution);
- Script 3 (plots for graphical visualisation of results).R: used for graphical visualisation of results, including (i) ordination plots (PCA, D-PGLS, pGLS, 2B-PLS), (ii) shape deformation 3D plots (PCA, D-PGLS, 2B-PLS), and (iii) trait-mapping on trees;
- Custom R functions.R: this script contains (i) custom R code used to retrieve results from D-PGLS objects, and (ii) custom functions from Motani and Schmitz (2011, Evolution) used to perform pFDA