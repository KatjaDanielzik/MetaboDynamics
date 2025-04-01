#@importFrom stats as.dendrogram 


# heatmap <- heatmap()
# 
# 
# # dendrogram
# plot(as.dendrogram(hclust))
# dendrogram <- recordPlot() # record plot for list
# ## or
# cluster::pltree()
# 
# 
# # PCA ?
# # base R
# prcomp() or princomp()
# PCA visualization with ggplot
# 
# 
# # plot dynamics
# temp <- temp%>%pivot_longer(cols=c("1","2","3","4"),names_to = "time.h",values_to = "mean_log_cpc_scaled")
# 
# lineplot <- ggplot(temp,aes(x=as.factor(as.numeric(time.h)),y=mean_log_cpc_scaled,group=metabolite))+
#   geom_line()+
#   theme_bw()+
#   facet_wrap(~cluster,scales="free_y")+
#   xlab("time (h)")+
#   ylim(c(-2.5,2.5))
# result <- list(heatmap = heatmap,
#                lineplot = lineplot)