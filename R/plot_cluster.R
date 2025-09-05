#' visualize clustering solution of cluster_dynamics()
#'
#'
#' @param data result of [cluster_dynamics()] function: either a list of dataframes or a SummarizedExperiment object
#'
#' @import ggplot2
#' @import tidyr
#' @import ggtree
#' @importFrom stats prcomp
#' @importFrom stats as.dendrogram
#' @importFrom stats order.dendrogram
#' @importFrom dendextend color_branches
#' @importFrom dendextend color_labels
#' @importFrom grDevices recordPlot
#' @importFrom graphics par
#' @importFrom graphics title
#'
#' @returns a list of plots. Per experimental condition: 1)a dendrogram of the
#' mean clustering solution, 2) a 'bubbletree': a phylogram with numbers on nodes
#' indicating in how many bootstraps of the posterior estimates the same clustering
#' solution was generated, 3) order of the tips in bubbletree: needed for matching lineplots and ORA,
#' 3) mean dynamics of clusters
#'
#' @export
#'
#' @seealso [cluster_dynamics()]
#'
#' @examples
#' data("longitudinalMetabolomics")
#' plot_cluster(longitudinalMetabolomics[, longitudinalMetabolomics$condition == "A"])

plot_cluster <- function(data){
  # Input checks
  if (!inherits(data, "list") && !inherits(data, "SummarizedExperiment")) {
    stop("'data' must be a list or a SummarizedExperiment object
         obtained by function cluster_dynamics()")
  }
  if (is(data, "list")) {
    lapply(data, function(data) {
      # check for dataframe format
      if (!inherits(data, "list")) {
        stop("'data' must be a list of lists obtained by cluster_dynamics() if it is not a SummarizedExperiment object")
      }
    })
  }

  # Data transformation
  if (is(data, "SummarizedExperiment")) {
    data <- metadata(data)[["cluster"]]
  }
  
  # color coded dendrograms of clustering solution of mean estimates
  dendrograms <- list()
  # plot dendrogram with clustering solution
  for (i in names(data[["cluster"]])){
      temp <- data[["cluster"]][[i]]
    # dendrogram
    dendro <- as.dendrogram(temp[["mean_dendro"]])
    # get color identifier based on cluster membership of metabolite
    colors_to_use <- as.numeric(as.data.frame(temp[["data"]][, c("metabolite", "cluster")])$cluster)
    # order by dendrogram
    colors_to_use <- colors_to_use[order.dendrogram(dendro)]
    dendro <- dendextend::color_branches(dendro, col = colors_to_use)
    dendro <- dendextend::color_labels(dendro, col = colors_to_use)
    # get labels
    labels <- temp[["mean_dendro"]]
    par(cex = 0.5)
    plot(dendro)
    par(cex = 1)
    title(main = paste0(
      i," dynamicTreeCut, method= ",
      labels$method
    ), ylab = paste0(labels$dist.method, " distance"))
    dendrograms[[i]] <- recordPlot()
  }
  
  # bubbletrees of bootstrapping
  trees <- list()
  for (i in names(data[["cluster"]])){
    temp <- data[["cluster"]][[i]]
    trees[[i]]<- ggtree(temp$mean_phylo,linetype='solid')+
      geom_point2(mapping = aes(subset=isTip==FALSE),size = 0.5, col = "black")+
      geom_tippoint(size = 2, fill = "white", shape = 21)+
      geom_tiplab(color='black', as_ylab = TRUE, align = TRUE)+
      layout_rectangular()+
      theme_bw(base_size = 10)+
      scale_x_continuous(labels = abs)+
      geom_nodelab(geom='text', color = "#4c4c4c" ,size = 2.75, hjust=-0.2,
                   mapping = aes(label=label,subset=isTip==FALSE))+
      ggtitle(paste0("condition ",i))
  }
  
  # plot dynamics as lineplots
  trees <- lapply(trees,revts)
  tips <- list()
  clusterplots <- list()
  lineplots <- list()
  for (i in names(trees)){
    tree <- trees[[i]]
    t  <- tree$data
    t <- t[order(t$y, decreasing = FALSE), ]
    tips[[i]] <- t$label[t$isTip==TRUE]
    temp <- data[["cluster"]][[i]]$data
    temp <- temp%>%pivot_longer(cols=-c(metabolite,KEGG,condition,cluster),names_to = "time",values_to = "mean")
    temp$metabolite <- factor(temp$metabolite,levels=tips[[i]])
    temp$cluster <- as.factor(temp$cluster)

    clusterplots[[i]]<-ggplot(temp,aes(y=metabolite,x=cluster,fill=cluster))+
        geom_tile()+
        scale_fill_viridis_d()+
        guides(col="cluster")+
        theme_bw()#+
        #scale_y_continuous(breaks=temp$tips_rank,labels = temp$metabolite)
  
    plots <- list()
    n_metabolites <- c()
    cluster_order <- unique(rev(temp[order(temp$metabolite),]$cluster))
    for (j in cluster_order){
      plots[[j]]<- temp%>%filter(cluster==j)%>%ggplot(aes(x=time,y=mean))+
      geom_line(aes(group=metabolite))+
      scale_color_viridis_d()+
      xlab("")+
      facet_grid(rows=vars(cluster))+
      theme_bw()+
      ylim(c(-2,2))+
      theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())
      summary <- temp%>%filter(cluster==j)%>%select(cluster,metabolite)%>%
        distinct()%>%
        group_by(cluster)%>%summarise(n_metabolite=n_distinct(metabolite))
      n_metabolites <- c(n_metabolites,summary$n_metabolite)
      plots[[paste0(j,"space")]] <- plot_spacer() # spacer plot to reduce space between plots
    }
    # assign heights according to number of metabolites
    heights <- as.vector(n_metabolites/sum(n_metabolites))*100 ## needs to be sorted!
    heights <- as.numeric(c(rbind(heights,-2.7))) # add spacing for spacer plots
    p <- Reduce("/",plots)
    lineplots[[i]] <- p + plot_layout(heights=heights)+plot_annotation("dynamics")+
      
  }
  
  i <- "A"
  trees[[i]]|clusterplots[[i]]|lineplots[[i]]
  
  return(list(dendrograms = dendrograms, trees = trees,
              clusterplots = clusterplots))
}
