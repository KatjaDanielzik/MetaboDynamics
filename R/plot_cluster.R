#' visualize clustering solution of cluster_dynamics()
#'
#'
#' @param data result of [cluster_dynamics()] function: either a list of data frames or a SummarizedExperiment object
#'
#' @import ggplot2
#' @import tidyr
#' @import ggtree
#' @import patchwork
#' @importFrom tidytree isTip
#'
#' @returns a list of plots. Per experimental condition: 1) a 'bubbletree': a phylogram with numbers on nodes
#' indicating in how many bootstraps of the posterior estimates the same clustering
#' solution was generated, 2) cluster affiliation of metabolites, 3) dynamics of metabolites per cluster,
#' 4) patchwork of 1-3, 5) order of the tips in bubbletree: needed for matching lineplots and ORA,
#' 3) mean dynamics of clusters
#'
#' @export
#'
#' @seealso [cluster_dynamics()]
#'
#' @examples
#' data("longitudinalMetabolomics")
#' plots <- plot_cluster(longitudinalMetabolomics[, longitudinalMetabolomics$condition == "A"])
#' 
#' plots[["trees"]][["A"]]

plot_cluster <- function(data){
  # Input checks
  if (!inherits(data, "list") && !inherits(data, "SummarizedExperiment")) {
    stop("'data' must be a list or a SummarizedExperiment object obtained by function cluster_dynamics()")
  }

  # Data transformation
  if (is(data, "SummarizedExperiment")) {
    data <- metadata(data)[["cluster"]]
  }
  if (!inherits(data, "list")){
    stop("clustering solution obtained by cluster_dynamics must be stored in metadata of SummarizedExperiment under 'cluster'")
  }
  # binding of global variables
  label <- NULL
  metabolite <- NULL
  condition <- NULL
  cluster <- NULL
  time <- NULL
  
  
  
  # bubbletrees of bootstrapping
  trees <- list()
  for (i in names(data)){
    temp <- data[[i]]
    trees[[i]]<- ggtree(temp$mean_phylo,linetype='solid')+
      geom_point2(mapping = aes(subset=isTip==FALSE),size = 0.5, col = "black")+
      geom_tippoint(size = 2, fill = "white", shape = 21)+
      geom_tiplab(color='black', as_ylab = TRUE, align = TRUE)+
      layout_rectangular()+
      theme_bw(base_size = 10)+
      scale_x_continuous(labels = abs)+
      geom_nodelab(geom='text', color = "#4c4c4c" ,size = 2.75, hjust=-0.2,
                   mapping = aes(label=label,subset=isTip==FALSE))
  }
  
  # plot dynamics as lineplots
  trees <- lapply(trees,revts)
  tips <- list()
  clusterplots <- list()
  lineplots <- list()
  cluster_order <- list()
  cluster_heights <- list()
  for (i in names(trees)){
    tree <- trees[[i]]
    t  <- tree$data
    t <- t[order(t$y, decreasing = FALSE), ]
    tips[[i]] <- t$label[t$isTip==TRUE]
    temp <- data[[i]]$data
    temp <- temp%>%pivot_longer(cols=-c(metabolite,condition,cluster),names_to = "time",values_to = "mean")
    temp$metabolite <- factor(temp$metabolite,levels=tips[[i]])
    temp$cluster <- as.factor(temp$cluster)

    clusterplots[[i]]<-ggplot(temp,aes(y=metabolite,x=cluster,fill=cluster))+
        geom_tile()+
        scale_fill_viridis_d()+
        guides(col="cluster")+
        ylab("")+
        theme_bw()
  
    plots <- list()
    n_metabolites <- c()
    cluster_order[[i]] <- unique(rev(temp[order(temp$metabolite),]$cluster)) # get factors 
    for (j in cluster_order[[i]]){
      plots[[j]]<- temp%>%filter(cluster==j)%>%ggplot(aes(x=time,y=mean))+
      geom_line(aes(group=metabolite))+
      scale_color_viridis_d()+
      xlab("")+
      facet_grid(rows=vars(cluster))+
      theme_bw()+
      ylim(c(-2,2))+
      theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())+
      geom_vline(aes(xintercept=as.factor(time)),col="grey",linetype="dashed")+
      geom_hline(aes(yintercept=0),col="grey",linetype="dashed")
      
      summary <- temp%>%filter(cluster==j)%>%select(cluster,metabolite)%>%
        distinct()%>%
        group_by(cluster)%>%summarise(n_metabolite=n_distinct(metabolite))
      n_metabolites <- c(n_metabolites,summary$n_metabolite)
      plots[[paste0(j,"space")]] <- plot_spacer() # spacer plot to reduce space between plots
    }
    # assign heights according to number of metabolites
    heights <- as.vector(n_metabolites/sum(n_metabolites))*100 ## needs to be sorted!
    heights <- as.numeric(c(rbind(heights,-2.7))) # add spacing for spacer plots
    cluster_heights[[i]] <- heights
    p <- Reduce("/",plots)
    lineplots[[i]] <- p + plot_layout(heights=heights)+plot_annotation("metabolite dynamics",
                                                                       "lines = metabolite, row panels = cluster,
    grey lines = time points")
  }
  
  patchwork <- list()
  for (i in names(data)){
  p <- trees[[i]]|clusterplots[[i]]|lineplots[[i]]
  patchwork[[i]] <- p+plot_annotation(paste0("condition ",i),
"plots: dendrogram with cluster probability, cluster affiliations of metabolites, metabolite dynamics in clusters")
  }
  
  return(list(trees = trees, clusterplots = clusterplots, lineplots=lineplots,
              patchwork = patchwork, cluster_order=cluster_order, cluster_heights = cluster_heights))
}
