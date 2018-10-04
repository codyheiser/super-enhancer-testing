# Data Generator and Testing Functions for SuperEnhancer R Package
#   define functions using code from AP to generate randomized B-globin expression data with different n
#   and iterations for testing of SuperEnhancer statistical model fitting package (Dukler, et al, 2017)
#
# C Heiser, 2018

require(reshape2)
require(plyr)
require(tidyverse)
require(ggpubr)
require(superEnhancerModelR)

# generate n random numbers that are Norm distributed with designated mean and sd
rnorm2 <- function(n,mean,sd,seed=NA){
  if(!is.na(seed)){
    set.seed(seed) # set seed for reproducible random number generation
  }
  return(mean+sd*scale(rnorm(n)))
} 

# generate data frame of random expression values using given means and sds
genBglobin <- function(n, wt.norm=F, set.seed = T, suppress.out=T){
  # n = number of replicates of expression data to randomly generate for each condition
  # wt.norm = generate random, normally-distributed points for WT set? default no.
  # suppress.out = write to .csv file? default no.
  
  if(set.seed){
    seed = 1:9
  }else{
    seed = rep(NA, 9)
  }
  
  ### From beta_globin-generateData.R, AP:
  ## creating data set for beta globin LCR based 
  ## on mean and standard deviation given in Bender 2001 and 2012
  if(wt.norm){
    WT <- rnorm2(n,1.00,0.05,seed[1]) # if specified, generate random WT data with mean of 1.0 and expected SD
  }else{
    WT <- rep(1,n) # otherwise, use default of 1.0 for all WT datapoints
  }
  D1 <- rnorm2(n,0.78,0.05,seed=seed[2]) # generated data for HS1 deletion 
  D2 <- rnorm2(n,0.59,0.04,seed=seed[3]) # generated data for HS2 deletion 
  D3 <- rnorm2(n,0.71,0.03,seed=seed[4]) # generated data for HS3 deletion 
  D4 <- rnorm2(n,0.81,0.05,seed=seed[5]) # generated data for HS4 deletion 
  D56 <- rnorm2(n,0.97,0.09,seed=seed[6]) # generated data for HS5-6 deletion 
  D14 <- rnorm2(n,0.60,0.03,seed=seed[7]) # generated data for HS1 and HS4 deletion 
  D12 <- rnorm2(n,0.39,0.01,seed=seed[8]) # generated data for HS1-2 deletion 
  D23 <- rnorm2(n,0.31,0.02,seed=seed[9]) # generated data for HS2-3 deletion 
  
  Bglobin <- data.frame(WT=WT, E1=D1, E2=D2, E3=D3, E4=D4, E56=D56, E12=D12, E14=D14, E23=D23)
  if(!suppress.out){
    write.csv(Bglobin, file = paste0('exprdata_',n,'reps_',gsub(Sys.time(),pattern = ' ',replacement = '_')), row.names = FALSE)
  }
  return(Bglobin)
}

# conditionally mutate rows of dataframe as part of dplyr::mutate function
mutate_cond <- function(.data, condition, ..., envir = parent.frame()){
  # This function performs a __mutate__ operation on a subset of rows of a dataframe without the need for __filter__ or __group_by__
  # and returns the output in place. Because of this, you cannot create new columns within this function as you can using normal
  # __dplyr::mutate__. Instead, ensure output columns are already initialized in the dataframe.
  condition <- eval(substitute(condition), .data, envir)
  .data[condition,] %>%
    mutate(...) -> .data[condition,]
  return(.data)
}

# put data from genBglobin into superEnhancerModelR format
reformatBglobin.superE <- function(df){
  df %>%
    melt(id.vars=c()) %>%
    dplyr::rename(condition=variable, expression=value) %>%
    mutate(E1=1,E2=1,E3=1,E4=1,E5=1,E6=1) %>%
    mutate_cond(str_detect(condition,'1'), E1=0) %>% 
    mutate_cond(str_detect(condition,'2'), E2=0) %>% 
    mutate_cond(str_detect(condition,'3'), E3=0) %>% 
    mutate_cond(str_detect(condition,'4'), E4=0) %>% 
    mutate_cond(str_detect(condition,'5'), E5=0) %>% 
    mutate_cond(str_detect(condition,'6'), E6=0) -> out
  return(out)
}

# generate model with given error and link functions
gen.model <- function(df, err, link, enhancer.formula = ~E1+E2+E3+E4+E5+E6, maxit=2000, 
                      actBounds=c(-150,150), errBounds=c(10^-3, 10^3), scaleBounds=c(10^-3, 10^3)){
  # df = data in superE format (i.e. genBglobin() %>% reformat_superE())
  # err = error function to use ('gaussian' or 'lognormal')
  # link = link function to use ('additive', 'exponential', 'logisitic')
  # enhancer.formula = interactions between variables to be modeled
  # maxit = maximum iterations of optimDE() to run
  # actBounds, errBounds, scaleBounds = options to pass to enhancerDataObject() that restrict search for coefficients
  
  start.time <- proc.time() # start timer
  
  if(maxit < 800){t <- 1}else{t <- 6} # determine parallelization based on number of iterations to run
  
  # create and optimize enhancer data model object
  enhance.obj <- enhancerDataObject(df[['expression']], df %>% select(-condition, -expression), formula(enhancer.formula), 
                                    errorModel = err, linkFunction = link, activityParameterBounds = actBounds, 
                                    errorParameterBounds = errBounds, scaleParameterBounds = scaleBounds) %>%
    optimDE(maxit=maxit, refine=T, threads=t, control=list(trace=round_any(maxit/10, 100, f = ceiling)))
  
  print(proc.time() - start.time) # report completion time
  
  return(enhance.obj)
}

# need below function for plotModel.CH
#   adapted from ndukler/superEnhancerModelR
errorIntervals.CH <- function(x,activity,quantiles){
  if(any(quantiles>0.5)){
    stop("Quantiles 0<quantiles<0.5")
  }
  err=list()
  express=predictExpression(x,activity)
  if(x@errorModel$type=="lognormal"){
    for(q in quantiles){
      err[[as.character(q)]]=data.frame(x=c(activity,rev(activity)),
                                        y=c(qlnorm(q, log(express), x@errorModel$value[1]),
                                            rev(qlnorm(1-q, log(express), x@errorModel$value[1]))),
                                        Quantile=paste0(q,"-",1-q))
    }
  } else if(x@errorModel$type=="gaussian"){
    for(q in quantiles){
      err[[as.character(q)]]=data.frame(x=c(activity,rev(activity)),
                                        y=c(qnorm(q, express, x@errorModel$value[1]),
                                            rev(qnorm(1-q, express, x@errorModel$value[1]))),
                                        Quantile=paste0(q,"-",1-q))
    }
  } else {
    stop("Unsupported error model")
  }
  return(do.call("rbind",err))
}

# preferred plotting options
#   call these by adding them (+) to a ggplot object
plot.opts <- list(
  theme_bw(),
  theme(legend.text=element_text(size=9, color = 'black'),
        axis.line = element_line(colour = 'black'),
        axis.title=element_text(size=12, color = 'black'),
        axis.text.x=element_text(size=10, color = 'black'),
        axis.text.y=element_text(size=10, color = 'black'),
        plot.title=element_text(size=12, color = 'black'))
)

# plot the output of the model
#   adapted from ndukler/superEnhancerModelR
plotModel.CH <- function(x){
  ## Get names of active enhancers for each experiment
  indiv.enhancers=grep(x = colnames(x@designMatrix),pattern = ")|:",invert = TRUE)
  rname=apply(x@designMatrix,1, function(z) paste(colnames(x@designMatrix)[indiv.enhancers][as.logical(z[indiv.enhancers])],collapse = "/"))
  rname[rname==""]="None"
  
  ## Get factor ordering of rows
  temp=data.frame(rname=rname,nenh=rowSums(x@designMatrix),stringsAsFactors = FALSE)
  
  ## Get activity values for actual oberservations
  real=data.frame(activity=computeActivity(x),observed=x@expressionData,enhancers=rname,stringsAsFactors = FALSE)
  real$enhancers=factor(real$enhancers,levels=unique(with(temp,temp[order(nenh,rname),])$rname))
  
  ## Compute expression values for intermediate activity to get smooth curve
  sim.act=seq(min(real$activity)-abs(min(real$activity))*0.1,max(real$activity)+abs(max(real$activity))*0.1,by=0.01)
  out=data.frame(activity=sim.act,expression=predictExpression(x,sim.act))
  
  err=errorIntervals.CH(x,sim.act,quantiles = 0.1)
  
  cc <- scales::seq_gradient_pal("light blue", "blue", "Lab")(seq(0,0.5,length.out=length(unique(err$Quantile))))
  
  g=ggplot2::ggplot()+
    ggplot2::geom_polygon(data=err,ggplot2::aes(x=x,y=y,fill=Quantile))+
    ggplot2::scale_fill_manual(values=cc)+
    ggplot2::geom_path(data=out,ggplot2::aes(activity,expression),color="black",size=1.5)+
    ggplot2::geom_point(data=real,ggplot2::aes(activity,observed,color=enhancers), size=2.5, alpha=0.6)+
    ggplot2::xlab("Activity/-Energy")+
    ggplot2::ylab("Expression")+
    labs(x="Activity/-Energy", y="Expression", color='Enhancers')+
    plot.opts
  
  return(g)
}

# function to test parameters of model
test.params <- function(df, errs, links, enhancer.formula = ~E1+E2+E3+E4+E5+E6, maxit=2000, 
                        actBounds=c(-150,150), errBounds=c(10^-3, 10^3), scaleBounds=c(10^-3, 10^3)){
  # df = data in superE format (i.e. genBglobin() %>% reformat_superE())
  # errs = error functions to use c('gaussian', 'lognormal')
  # links = link functions to use c('additive', 'exponential', 'logisitic')
  # enhancer.formula = interactions between variables to be modeled
  # maxit = maximum iterations of optimDE() to run
  # actBounds, errBounds, scaleBounds = options to pass to enhancerDataObject() that restrict search for coefficients
  
  start.time <- proc.time() # start timer
  
  bic <- NULL # initiate object to track BICs
  mods <- list() # initiate empty list for dumping models into
  plts <- list() # initiate empty list for dumping plots into
  
  for(link.function in links){
    for(err.function in errs){
      
      mod <- NULL
      try(mod <- gen.model(df, err = err.function, link = link.function, enhancer.formula = enhancer.formula, maxit = maxit, 
                           actBounds = actBounds, errBounds = errBounds, scaleBounds = scaleBounds)) # model with given parameters
      if(is.null(mod)){next} # if model fails, move on to next link-error combination
      
      mods[[length(mods) + 1]] <- mod # add model to list
      plts[[length(plts) + 1]] <- plotModel.CH(mod)+labs(title = paste0(link.function,'/',err.function)) # add plot of results
      
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
  
  print(proc.time()-start.time)
  return(list(bic,mods,plts)) # return as list of objects
}

# define function to plot multiple complete plot objects on one image
sum.fig.superE <- function(plotlist, bic.vals = NULL){
  # plotlist = list of ggplot objects generated from superEnhancerModelR with all possible link/error function combinations 
  #   (see superE_datagen_CH.R)
  # bic.vals = table of BIC values from test.params(), optional
  
  # clean plots for arranging in figure
  clean.plts <- lapply(plotlist, FUN = function(x){return(x+labs(x=NULL,y=NULL,color='Enhancers')+theme_pubr())})
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

# extract link function coefficient values and generate plots

