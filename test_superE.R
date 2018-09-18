# Testing SuperEnhancer R Package
#   for testing of SuperEnhancer statistical model fitting package (Dukler, et al, 2017)
#
# C Heiser, 2018

rm(list=ls()) # clear workspace
suppressPackageStartupMessages(require(argparse))
suppressPackageStartupMessages(source('utilityfunctions_superE.R')) # source functions needed to perform tests

# create parser object
parser <- ArgumentParser()
# import options
parser$add_argument('data',
                    help='Path to data in superEnhancerModelR format as .csv file.')
parser$add_argument('outputlocation', 
                    help='Path to directory to save outputs.')
parser$add_argument('-i', '--iterations', nargs='+', 
                    help='Number of iterations to perform in optimDE(). Can be list of values.')
parser$add_argument('-ab', '--activitybounds',  default='c(10^-3, 10^3)',
                    help='Activity bounds to pass to superEnhancerDataObject() function.')
parser$add_argument('-eb', '--errorbounds',  default='c(10^-3, 10^3)',
                    help='Error parameter bounds to pass to superEnhancerDataObject() function.')
parser$add_argument('-sb', '--scalebounds',  default='c(10^-3, 10^3)',
                    help='Scale parameter bounds to pass to superEnhancerDataObject() function.')
# get command line options, if help encountered print help and exit,
#   otherwise if options not found on command line, set defaults
args <- parser$parse_args()
# read data into df
datain <- read.csv(args$data)

# define parameters to test
error.models <- c('gaussian','lognormal')
link.functions <- c('additive','exponential','logistic')

# define testing function
test.model <- function(df, optim.iter, out = 'outputs/', ...){
  # df = data in superE format
  # optim.iter = total number of iteration for optimization function to perform. can be list. 
  # out = path to output directory
  # ... = additional parameters to pass to enhancerDataObject()
  
  master.out <- data.frame() # initiate df for dumping fit parameters and BICs into
  
  for (iterations in optim.iter) {
      # print some useful stuff to console
      print(paste0('Performing ', iterations, ' iterations on ', args$data))
    
      # run test using all six link/error function combos
      result <- test.params(df, error.models, link.functions, maxit = iterations, ...) # build models using all error/link function combinations
        
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
        temp.link <- cbind(temp.link, data.frame(optim.iter = iterations))
        # add linkFunction to master frame
        master.link <- rbind.fill(master.link, temp.link)
      }
      # merge link df with BIC df and append to master.out
      master.out <- rbind.fill(master.out, merge(master.link, result[[1]], by = c('link','error'))) 
        
      # generate pretty figure and save to .pdf file
      figure <- sum.fig.superE(result[[3]], bic.vals = result[[1]])
      ggsave(figure, filename = paste0(out, 'superE_Summary_', iterations, 'iter_', Sys.Date(),'.pdf'), 
             device = 'pdf', width = 12, height = 8, units = 'in')
  }
  write.csv(master.out, file = paste0(out,'superE_Results_',Sys.Date(),'.csv'), row.names = F)
}

test.model(df = datain, optim.iter = as.numeric(args$iterations), 
           out = args$outputlocation, activityParameterBounds = eval(parse(text=args$activitybounds)), 
           errorParameterBounds = eval(parse(text=args$errorbounds)), scaleParameterBounds = eval(parse(text=args$scalebounds)))