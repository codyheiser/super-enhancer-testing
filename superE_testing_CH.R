# Testing SuperEnhancer R Package
#   use code from AP to generate randomized B-globin expression data with different n and iterations 
#   for testing of SuperEnhancer statistical model fitting package (Dukler, et al, 2018)
#
# C Heiser, 2018

rm(list=ls()) # clear workspace

require(plyr)
require(testthat)
require(ggpubr)
source('superE_datagen_CH.R') # source functions and stuff needed to perform tests

# define function to plot multiple complete plot objects on one image
sum.fig.superE <- function(plotlist, bic.vals = NULL){
  # plotlist = list of ggplot objects generated from superEnhancerModelR with all possible link/error function combinations 
  #   (see superE_datagen_CH.R)
  # bic.vals = table of BIC values from test.params(), optional
  
  # clean plots for arranging in figure
  clean.plts <- lapply(plotlist, FUN = function(x){return(x+labs(title=NULL,x=NULL,y=NULL,color='Enhancers')+theme_pubr())})
  # arrange model plots into figure with common legend and clean graphs
  fig <- ggarrange(plotlist = clean.plts, ncol = 2, nrow = 3, common.legend = T, legend = 'right') %>%
    annotate_figure(left = text_grob('Expression', rot = 90, size = 14), bottom = text_grob('Activity/-Energy', size = 14))
  # include bic plot if available
  if(!is.null(bic.vals)){
    bic.plt <- ggplot(data=bic.vals, aes(x=link, y=bic, fill=error))+
      geom_bar(stat='identity', position = position_dodge(width = 1))+
      labs(x=NULL, y='BIC', fill='Error Function')+
      plot.opts
    rel.bic.plt <- ggplot(data=bic.vals, aes(x=link, y=rel.bic, fill=error))+
      geom_bar(stat='identity', position = position_dodge(width = 1))+
      labs(x=NULL, y='Relative BIC', fill='Error Function')+
      plot.opts
    # add bic to final figure
    fig <- ggarrange(fig, ggarrange(bic.plt, rel.bic.plt, NULL, nrow = 3, labels = c('B','C')), ncol = 2, widths = c(2,1), labels = 'A')
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
ggsave(sum.fig.wap, filename = 'outputs/wap_summary.pdf', device = 'pdf', width = 12, height = 8, units = 'in')


#############################################################################################################################################
# Example alpha-globin data from super-enhancer-code and Dukler, et al (2017)
Aglobin <- read.csv('inputs/alpha_globin_data.csv') # load alpha-globin data
# TODO: reshape data

# generate results for all link/error function combinations
Aglobin_results <- test.params(Aglobin, error.models, link.functions, enhancer.formula = ~E1+E2+E3) # TODO: get proper enhancer formula

# generate figure
sum.fig.Aglobin <- sum.fig.superE(Aglobin_results[[3]], bic.vals = Aglobin_results[[1]])
ggsave(sum.fig.Aglobin, filename = 'outputs/Aglobin_summary.pdf', device = 'pdf', width = 12, height = 8, units = 'in')


#############################################################################################################################################
# do some testing with contrived  Bglobin data:
test.Bglobin <- function(expr.reps, wt.norm, optim.iter, out = 'outputs/', ...){
  # expr.reps = number of replicates of expression data for each enhancer condition. can be list. 
  # wt.norm = if TRUE, generate normally-distributed datapoints around 1 to represent WT expression. can be list. 
  # optim.iter = total number of iteration for optimization function to perform. can be list. 
  # out = path to output directory
  # ... = additional parameters to pass to enhancerDataObject()
  
  bglobin.results <- list()
  bglobin.figures <- list()
  master.out <- data.frame() # initiate df for dumping fit parameters and BICs into
  
  for (norm.strategy in wt.norm) {
    for (reps in expr.reps) {
      for (iterations in optim.iter) {
        # print some useful stuff to console
        print(paste0('Performing ', iterations, ' iterations on ', reps, ' expression data points with WT normalization == ', norm.strategy))
        
        # generate Bglobin data and put into superE format
        bglobin <- genBglobin(reps, norm.strategy) %>% reformatBglobin.superE()
        
        # run test on Bglobin data using all six link/error function combos
        result <- test.params(bglobin, error.models, link.functions, maxit = iterations, ...) # build models using all error/link function combinations
        
        # append important info to master.out df
        master.link <- data.frame()
        for (model in 1:length(result[[2]])) {
          # get linkFunction$value into df format
          temp.link <- as.data.frame(result[[2]][[model]]@linkFunction$value) %>%
            mutate(variable = rownames(as.data.frame(result[[2]][[model]]@linkFunction$value))) %>%
            dplyr::rename(value = 1) %>%
            spread(variable, value) %>%
            # append link and error function info
            mutate(link = result[[2]][[model]]@linkFunction$type, 
                   error = result[[2]][[model]]@errorModel$type, 
                   error.sd = result[[2]][[model]]@errorModel$value[['sd']])
          temp.link <- cbind(temp.link, data.frame(norm.strategy = norm.strategy, expr.reps = reps, optim.iter = iterations))
          # add linkFunction to master frame
          master.link <- rbind.fill(master.link, temp.link)
        }
        # merge link df with BIC df and append to master.out
        master.out <- rbind.fill(master.out, merge(master.link, result[[1]], by = c('link','error'))) 
        
        # generate pretty figure and save to .pdf file
        figure <- sum.fig.superE(result[[3]], bic.vals = result[[1]])
        ggsave(figure, filename = paste0(out, 'Bglobin_', reps, '_', norm.strategy, '_', iterations, '.pdf'), 
               device = 'pdf', width = 12, height = 8, units = 'in')
        
        # append to lists for massive output
        bglobin.results[[length(bglobin.results) + 1]] <- result # add model result to list
        bglobin.figures[[length(bglobin.figures) + 1]] <- figure # add model figure to list
      }
    }
  }
  write.csv(master.out, file = paste0(out,'results_',Sys.Date(),'.csv'), row.names = F)
  return(list(master.out, bglobin.results, bglobin.figures))
}

master2 <- test.Bglobin(expr.reps = 10, wt.norm = FALSE, optim.iter = c(500, 2000, 4000, 10000), activityParameterBounds = c(10^-2, 10^2),
                        errorParameterBounds = c(10^-2, 10^2), scaleParameterBounds = c(10^-2, 10^2), 
                        out = '~/Dropbox/_Venters_Lab_Resources/3_Rotation_Students/4_Cody/superE/parameter_testing_091718/')
