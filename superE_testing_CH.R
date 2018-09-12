# Testing SuperEnhancer R Package
#   use code from AP to generate randomized B-globin expression data with different n and iterations 
#   for testing of SuperEnhancer statistical model fitting package (Dukler, et al, 2018)
#
# C Heiser, 2018

source('superE_datagen_CH.R') # source functions and stuff needed to perform tests
set.seed(18) # set seed for reproducible random number generation

# define function to plot multiple complete plot objects on one image
sum.fig.superE <- function(plotlist, bic.vals = NULL){
  # plotlist = list of ggplot objects generated from superEnhancerModelR with all possible link/error function combinations 
  #   (see superE_datagen_CH.R)
  # bic.plot = table of BIC values from test.params(), optional
  
  # clean plots for arranging in figure
  clean.plts <- lapply(plotlist, FUN = function(x){return(x+labs(title=NULL,x=NULL,y=NULL,color='Enhancers')+theme_pubr())})
  # arrange model plots into figure with common legend and clean graphs
  fig <- ggarrange(plotlist = list(clean.plts[[1]], clean.plts[[4]], clean.plts[[2]], clean.plts[[5]], clean.plts[[3]], 
                                   clean.plts[[6]]), ncol = 2, nrow = 3, common.legend = T, legend = 'right') %>%
    annotate_figure(left = text_grob('Expression', rot = 90, size = 14), bottom = text_grob('Activity/-Energy', size = 14))
  # include bic plot if available
  if(!is.null(bic.vals)){
    bic.plt <- ggplot(data=bic.vals, aes(x=link, y=bic, fill=error))+
      geom_bar(stat='identity', position = position_dodge(width = 1))+
      labs(x=NULL, y='BIC', fill='Error Function', title='BIC for Superenhancer Models')+
      plot.opts
    # add bic to final figure
    fig <- ggarrange(fig, ggarrange(bic.plt, NULL, NULL, nrow = 3), ncol = 2, widths = c(2,1), labels = 'AUTO')
  }
  return(fig)
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

# generate figure
sum.fig.wap <- sum.fig.superE(wap_results[[3]], bic.vals = wap_results[[1]])
ggsave(sum.fig.wap, filename = 'wap_summary.pdf', device = 'pdf', width = 12, height = 8, units = 'in')


#############################################################################################################################################
# 10 reps, WT fixed
fix_10 <- genBglobin(10, wt.norm = F) %>% reformat_superE()
fix_10_results <- test.params(fix_10, error.models, link.functions)

# generate figure
sum.fig.fix_10 <- sum.fig.superE(fix_10_results[[3]], bic.vals = fix_10_results[[1]])
ggsave(sum.fig.fix_10, filename = 'Bglobin_10_fix_summary.pdf', device = 'pdf', width = 12, height = 8, units = 'in')

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
# 100 reps, WT fixed
fix_100 <- genBglobin(100, wt.norm = F) %>% reformat_superE()
fix_100_results <- test.params(fix_100, error.models, link.functions)

# generate figure
sum.fig.fix_100 <- sum.fig.superE(fix_100_results[[3]], bic.vals = fix_100_results[[1]])
ggsave(sum.fig.fix_100, filename = 'Bglobin_100_fix_summary.pdf', device = 'pdf', width = 12, height = 8, units = 'in')

# save plots for each model individually
i <- 1
for(plt in fix_100_results[[2]]){
  ggsave(plt, filename = paste0('wtfixed_100reps/Bglobin_fix100_',i,'.svg'), device = 'svg')
  png(file = paste0('wtfixed_100reps/Bglobin_fix100_',i,'.png'), width = 5, height = 3.6, units = 'in', res = 500)
  plot(plt)
  dev.off()
  i <- i + 1
}


#############################################################################################################################################
# 10 reps, WT norm
norm_10 <- genBglobin(10, wt.norm = T) %>% reformat_superE()
norm_10_results <- test.params(norm_10, error.models, link.functions)

# generate figure
sum.fig.norm_10 <- sum.fig.superE(norm_10_results[[3]], bic.vals = norm_10_results[[1]])
ggsave(sum.fig.norm_10, filename = 'Bglobin_10_norm_summary.pdf', device = 'pdf', width = 12, height = 8, units = 'in')

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
# 100 reps, WT norm
norm_100 <- genBglobin(100, wt.norm = T) %>% reformat_superE()
norm_100_results <- test.params(norm_100, error.models, link.functions)

# generate figure
sum.fig.norm_100 <- sum.fig.superE(norm_100_results[[3]], bic.vals = norm_100_results[[1]])
ggsave(sum.fig.norm_100, filename = 'Bglobin_100_norm_summary.pdf', device = 'pdf', width = 12, height = 8, units = 'in')

# save plots for each model individually
i <- 1
for(plt in norm_100_results[[2]]){
  ggsave(plt, filename = paste0('wtnorm_100reps/Bglobin_norm100_',i,'.svg'), device = 'svg')
  png(file = paste0('wtnorm_100reps/Bglobin_norm100_',i,'.png'), width = 5, height = 3.6, units = 'in', res = 500)
  plot(plt)
  dev.off()
  i <- i + 1
}

