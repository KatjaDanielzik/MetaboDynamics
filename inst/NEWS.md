# MetaboDynamics 1.1.7
- altered plot_cluster and plot_ORA visualization for patchwork option: easier visualization

# MetaboDynamics 1.1.6
- get_ORA_annotations function (retrieving KEGG IDs hierarchies) added back to package
- no more samples in estimates dynamics as probability of clustering solution will be implemented with a bubbletree in cluster_dynamics
- estimates now returns a list with: estimated metabolite abundance (mu), 
  estimated standard deviation of metabolite abundance (sigma), 
  estimated pooled standard deviation per metabolite and dose (lambda),
  differences in metabolite abundances between time points, 
  euclidean distance between metabolite dynamics of different conditions
- plot estimates additionally visualized the euclidean distance between conditions of metabolite specific dynamics
- cluster_dynamics provides clustering solution of mean estimates of mu as well as bootstraps clustering solutions
- plot cluster provides bubbletree, cluster identity, dynamics plots as well as patchwork plot
- plot_ORA has now option to be added to bubbletree obtained by plot_cluster
- all functions require named columns "metabolite","condition","time", "KEGG" (for ORA))
- ORA_hypergeometric() requires now a data frame annotationg metabolites to KEGG IDs

# MetaboDynamics 1.1.5
- internal adding of ANOVA model with euclidean distance estimation between doses and ANOVA model that integrates cell count uncertainty
- internal adding of simulated cell counts to data("longitudinalMetabolomics")

# MetaboDynamics 1.1.4
- differences between time points are now ordered and return is a list of plots 

# MetaboDynamics 1.1.2
contains a vignette describing the package workflow with a data frame input

# MetaboDynamics 1.0.2
Minor bug fix in vignette that caused errors in package checks

# MetaboDynamics 1.0.1
bug fix: fit_dynamics_model() can now correctly handle provided variable names

# MetaboDynamics 1.0.0
Bioconductor Release 3_21

# MetaboDynamics 0.99.26
* deprecate function get_ORA_annotations, needed data for ORA_hypergeometric
is included in package ("modules_compounds")

# MetaboDynamics 0.99.25
* comparison of clusters including visualization now handles number of clusters >10 correctly
* altered vignette including new clustering functions

# MetaboDynamics 0.99.24

* Bug fix: estimates_dynamics handles now differences between timepoints correctly
* added functions: cluster_dynamics: convenient clustering function wrapper using dynamictreecut package
                    plot_cluster: cluster visualization

# MetaboDynamics 0.99.22

* Bug fix: fit_dynamics allows number of time points that are not equal to 4

# MetaboDynamics 0.99.20

* Bug fix: can now handle tibbles as input (along with data frames and summarizedExperiments)

# Metabodynamics 0.99.18

* Accepted to Bioconductor https://www.bioconductor.org/packages/devel/bioc/html/MetaboDynamics.html 

# MetaboDynamics 0.99.0

* Initial Bioconductor submission.


