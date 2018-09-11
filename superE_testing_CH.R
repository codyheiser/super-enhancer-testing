# Testing SuperEnhancer R Package
#   use code from AP to generate randomized B-globin expression data with different n and iterations 
#   for testing of SuperEnhancer statistical model fitting package (Dukler, et al, 2018)
#
# C Heiser, 2018

source('superE_datagen_CH.R') # source functions and stuff needed to perform tests
set.seed(18) # set seed for reproducible random number generation

# define function to plot multiple complete plot objects on one image
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL){
  # Multiple plot function
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if(is.null(layout)){
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if(numPlots==1){
    print(plots[[1]])
  }else{
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for(i in 1:numPlots){
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }}
}

# define parameters to test
error.models <- c('gaussian','lognormal')
link.functions <- c('additive','exponential','logistic')


#############################################################################################################################################
# Example wap data from superEnhancerModelR package
data('wap') # load wap data
wap.design = wap[,c(3:5)]

# generate results for all link/error function combinations
wap_results <- test.params(wap, error.models, link.functions, enhancer.formula = ~E1+E2+E3)

# plot BIC values
wap.bic.plt <- ggplot(data=wap_results[[1]], aes(x=link, y=bic, fill=error))+
  geom_bar(stat='identity', position = position_dodge(width = 1))+
  labs(x=NULL, y='BIC', fill='Error Function', title='BIC for Superenhancer Models', subtitle='example wap data')+
  plot.opts
ggsave(wap.plt, filename = 'BIC_wap.pdf', width = 6, height = 4, units = 'in')

# plot everything at once
wap.sum.plt <- multiplot(plotlist = lapply(wap_results[[2]], FUN = function(x){
  return(x+theme(legend.position = 'none', plot.title = element_text(size=10), axis.title = element_blank()))
}), cols = 6)
ggsave(wap.sum.plt, filename = 'summary_wap.pdf', width = 18, height = 6, units = 'in')


#############################################################################################################################################
# 10 reps, WT norm
norm_10 <- genBglobin(10, wt.norm = T) %>% reformat_superE()
norm_10_results <- test.params(norm_10, error.models, link.functions)
# plot BIC values
png(file = 'wtnorm_10reps/BIC_Bglobin_norm_10.png', width = 600, height = 300, res = 100)
ggplot(data=norm_10_results[[1]], aes(x=link, y=bic, fill=error))+
  geom_bar(stat='identity', position = position_dodge(width = 1))+
  labs(x=NULL, y='BIC', fill='Error Function', title='BIC for Superenhancer Models', subtitle='10 replicates, normally-distributed WT data')+
  plot.opts
dev.off()
# plot everything at once
pdf(file = 'wtnorm_10reps/summary_Bglobin_norm_10.pdf', width = 18, height = 6)
multiplot(plotlist = lapply(norm_10_results[[2]], FUN = function(x){
  return(x+theme(legend.position = 'none', plot.title = element_text(size=10), axis.title = element_blank()))
}), cols = 6)
dev.off()
# save plots for each model individually
i <- 1
for(plt in norm_10_results[[2]]){
  ggsave(plt, filename = paste0('wtnorm_10reps/Bglobin_norm10_',i,'.svg'), device = 'svg')
  png(file = paste0('wtnorm_10reps/Bglobin_norm10_',i,'.png'), width = 5, height = 3.6, units = 'in', res = 500)
  plot(plt)
  dev.off()
  i <- i + 1
}

#############################################################################################################################################
# 10 reps, WT fixed
fix_10 <- genBglobin(10, wt.norm = F) %>% reformat_superE()
fix_10_results <- test.params(fix_10, error.models, link.functions)
# plot BIC values
png(file = 'wtfixed_10reps/BIC_Bglobin_fix_10.png', width = 1200, height = 800, res = 150)
ggplot(data=fix_10_results[[1]], aes(x=link, y=bic, fill=error))+
  geom_bar(stat='identity', position = position_dodge(width = 1))+
  labs(x=NULL, y='BIC', fill='Error Function', title='BIC for Superenhancer Models', subtitle='10 replicates, fixed WT data')+
  plot.opts
dev.off()
# plot everything at once
pdf(file = 'wtfixed_10reps/summary_Bglobin_fix_10.pdf', width = 18, height = 6)
multiplot(plotlist = lapply(fix_10_results[[2]], FUN = function(x){
  return(x+theme(legend.position = 'none', plot.title = element_text(size=10), axis.title = element_blank()))
}), cols = 6)
dev.off()
# save plots for each model individually
i <- 1
for(plt in fix_10_results[[2]]){
  ggsave(plt, filename = paste0('wtfixed_10reps/Bglobin_fix10_',i,'.svg'), device = 'svg')
  png(filename = paste0('wtfixed_10reps/Bglobin_fix10_',i,'.png'), width = 5, height = 3.6, units = 'in', res = 500)
  plot(plt)
  dev.off()
  i <- i + 1
}

#############################################################################################################################################
# 100 reps, WT norm
# Error in optim(par = c(-2.41105048323704, 0.0545654429076843, 0.81681304851154,  : 
#                          L-BFGS-B needs finite values of 'fn'
norm_100 <- genBglobin(100, wt.norm = T) %>% reformat_superE()
norm_100_results <- test.params(norm_100, error.models, link.functions)
# plot BIC values
png(file = 'wtnorm_100reps/BIC_Bglobin_norm_100.png', width = 1200, height = 800, res = 150)
ggplot(data=norm_100_results[[1]], aes(x=link, y=bic, fill=error))+
  geom_bar(stat='identity', position = position_dodge(width = 1))+
  labs(x=NULL, y='BIC', fill='Error Function', title='BIC for Superenhancer Models', subtitle='100 replicates, normally-distributed WT data')+
  plot.opts
dev.off()
# plot everything at once
pdf(file = 'wtnorm_100reps/summary_Bglobin_norm_100.pdf', width = 18, height = 6)
multiplot(plotlist = lapply(norm_100_results[[2]], FUN = function(x){
  return(x+theme(legend.position = 'none', plot.title = element_text(size=10), axis.title = element_blank()))
}), cols = 6)
dev.off()
# save plots for each model individually
i <- 1
for(plt in norm_100_results[[2]]){
  ggsave(plt, filename = paste0('wtnorm_100reps/Bglobin_norm100_',i,'.svg'), device = 'svg')
  png(file = paste0('wtnorm_100reps/Bglobin_norm100_',i,'.png'), width = 5, height = 3.6, units = 'in', res = 500)
  plot(plt)
  dev.off()
  i <- i + 1
}

#############################################################################################################################################
# 100 reps, WT fixed
# TODO: exponential with lognormal err function (05Sep18)
# Error in optim(par = c(-2.41105048323704, 0.0545654429076843, 0.81681304851154,  : 
#                          L-BFGS-B needs finite values of 'fn'
fix_100 <- genBglobin(100, wt.norm = F) %>% reformat_superE()
fix_100_results <- test.params(fix_100, error.models, link.functions)
# plot BIC values
png(file = 'wtfixed_100reps/BIC_Bglobin_fix_100.pdf', width = 1200, height = 800, res = 150)
ggplot(data=fix_100_results[[1]], aes(x=link, y=bic, fill=error))+
  geom_bar(stat='identity', position = position_dodge(width = 1))+
  labs(x=NULL, y='BIC', fill='Error Function', title='BIC for Superenhancer Models', subtitle='100 replicates, fixed WT data')+
  plot.opts
dev.off()
# plot everything at once
pdf(file = 'wtfixed_100reps/summary_Bglobin_fix_100.pdf', width = 18, height = 6)
multiplot(plotlist = lapply(fix_100_results[[2]], FUN = function(x){
  return(x+theme(legend.position = 'none', plot.title = element_text(size=10), axis.title = element_blank()))
}), cols = 6)
dev.off()
# save plots for each model individually
i <- 1
for(plt in fix_100_results[[2]]){
  ggsave(plt, filename = paste0('wtfixed_100reps/Bglobin_fix100_',i,'.svg'), device = 'svg')
  png(file = paste0('wtfixed_100reps/Bglobin_fix100_',i,'.png'), width = 5, height = 3.6, units = 'in', res = 500)
  plot(plt)
  dev.off()
  i <- i + 1
}
