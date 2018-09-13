# Super Enhancer Model Testing
#   define functions to test super enhancer statistical model fitting package (Dukler, et al. 2017)
#
# C Heiser, 2018

rm(list=ls()) # clear workspace

require(reshape2)
require(dplyr)
require(stringr)
require(tidyr)
require(ggplot2)
require(ggpubr)
require(superEnhancerModelR)

#############################################################################################################################################
# define functions

# preferred plotting options
# call these by adding them (+) to a ggplot object
plot.opts <- list(
  theme_bw(),
  theme(legend.text=element_text(size=9, color = 'black'),
        axis.line = element_line(colour = 'black'),
        axis.title=element_text(size=12, color = 'black'),
        axis.text.x=element_text(size=10, color = 'black'),
        axis.text.y=element_text(size=10, color = 'black'),
        plot.title=element_text(size=12, color = 'black'))
)

# generate model with given error and link functions
gen.model <- function(df, err, link, enhancer.formula = ~E1+E2+E3+E4+E5+E6, ...){
  # df = data in superE format 
  # err = error function to use ('gaussian' or 'lognormal')
  # link = link function to use ('additive', 'exponential', 'logisitic')
  # enhancer.formula = interactions between variables to be modeled
  # ... = options to pass to optimDE
  
  start.time <- proc.time() # start timer
  
  expr <- df[['expression']] # pull out vector of expression data
  design <- df %>% select(-condition, -expression) # pull out design matrix
  actFun <- formula(enhancer.formula) # create activity function
  enhance.obj <- enhancerDataObject(expr, design, actFun, errorModel = err, linkFunction = link) %>%
    optimDE(maxit=2000, refine=T, threads=6, control=list(trace=500), ...)
  
  print(proc.time() - start.time) # report completion time
  return(enhance.obj)
}

# function to test parameters of model
test.params <- function(df, errs, links, ...){
  # df = data in superE format (i.e. genBglobin() %>% reformat_superE())
  # errs = error functions to use c('gaussian', 'lognormal')
  # links = link functions to use c('additive', 'exponential', 'logisitic')
  
  bic <- NULL # initiate object to track BICs
  mods <- list() # initiate empty list for dumping models into
  plts <- list() # initiate empty list for dumping plots into
  resids <- list() # initiate empty list for dumping residual plots into
  
  for(link.function in links){
    for(err.function in errs){
      
      mod <- NULL
      try(mod <- gen.model(df, err.function, link.function, ...)) # model with given parameters
      if(is.null(mod)){next} # if model fails, move on to next link-error combination
      
      mods[[length(mods) + 1]] <- mod # add model to list
      plts[[length(plts) + 1]] <- plotModel(mod)+labs(title = paste0(link.function,'/',err.function)) # add plot of results
      resids[[length(resids) + 1]] <- plotResiduals(mod)+labs(title = paste0(link.function,'/',err.function,' residuals'))+plot.opts # add plot of residuals
      
      print(paste0(link.function,' / ',err.function,' BIC: ',bic(mod))) # print BIC
      
      # append BIC information to df for output
      if(is.null(bic)){
        bic <- data.frame(link = link.function, error = err.function, bic = bic(mod))
      }else{
        bic <- rbind(bic, data.frame(link = link.function, error = err.function, bic = bic(mod)))
      }
    }
  }
  bic$rel.bic <- bic$bic - (bic %>% filter(link=='additive', error=='gaussian'))$bic # calculate BICs relative to additive/gaussian combo
  return(list(bic,mods,plts,resids)) # return as list of objects
}

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


#############################################################################################################################################
# test example wap data from superEnhancerModelR package

# define parameters to test
error.models <- c('gaussian','lognormal')
link.functions <- c('additive','exponential','logistic')

data('wap') # load wap data
wap.design = wap[,c(3:5)]

# generate results for all link/error function combinations
wap_results <- test.params(wap, error.models, link.functions, enhancer.formula = ~E1+E2+E3)

# generate figure
sum.fig.wap <- sum.fig.superE(wap_results[[3]], bic.vals = wap_results[[1]])
ggsave(sum.fig.wap, filename = 'wap_summary.pdf', device = 'pdf', width = 12, height = 8, units = 'in')
