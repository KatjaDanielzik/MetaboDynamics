# MetaboDynamics 2.1.104
- bug fix in estimates_dynamics leading to mixed up condition annotations to
model fit results

# MetaboDynamics 2.1.103
- bug fix in estimates_dynamics leading to mixed up metabolite annotations to 
model fit results

# MetaboDynamics 2.1.101
- bug fix in fit_dynamics_model()
- new message in diagnostics_dynamics

# MetaboDynamics 2.1.100
- bug fix in plot_PPC for raw_plus_counts_model
- estimates_dynamics does not return 'sigma' or 'lambda' anymore

# MetaboDynamics 2.1.99
- added model_option in fit_dynamics_model(): either "sd_per_time_point" (only option before)
or "sd_per_condition" (less sigmas, more robust with few replicates)

# MetaboDynamics 2.1.2
- bug fix in plot_cluster: correct sorting of time points

# MetaboDynamics 2.1.1
- cluster_dynamics bug fix: can handle all metabolite names

# MetaboDynamics 1.99.12
- contains a **vignette describing the package workflow with a data frame input**
- altered plot_cluster and plot_ORA visualization for patchwork option
enabling easier visualization
- **all functions** require columns names "metabolite","condition","time" for consistency throughout the package
- **get_ORA_annotations()** function (retrieving KEGG IDs hierarchies) added back to package
- no more samples in estimates dynamics as probability of clustering solution is 
implemented now as bootstrapping from model posterior in function cluster_dynamics
- **estimates_dynamics()** now returns a list with: estimated metabolite abundance (mu), 
  estimated standard deviation of metabolite abundance (sigma), 
  estimated pooled standard deviation per metabolite and dose (lambda),
  differences in metabolite abundances between time points, 
  euclidean distance between metabolite dynamics of different conditions
- **plot_estimates()** additionally visualizes the euclidean distance between conditions of metabolite specific dynamics
- **cluster_dynamics()** provides clustering solution of mean estimates of mu as well as bootstraps clustering solutions
- **plot_cluster** provides bubbletree, cluster identity, dynamics plots as well as patchwork plot (combining bubbletree, cluster identity and dynamics plots)
- **plot_ORA()** has now option to be added to bubbletree obtained by plot_cluster
- **ORA_hypergeometric()** requires now a data frame annotationg metabolites to KEGG IDs (column names "metabolite" and "KEGG")
- internal adding of ANOVA model with euclidean distance estimation between doses and ANOVA model that integrates cell count uncertainty
- internal adding of simulated cell counts to data("longitudinalMetabolomics")
- differences between time points are now ordered and return is a list of plots 

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


