# Testing SuperEnhancer R Package
#   use code from AP to generate randomized B-globin expression data with different n and iterations 
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
parser$add_argument('-i', '--iterations', nargs='+', default='2000',
                    help='Number of iterations to perform in optimDE(). Can be list of values. Default 2000.')
parser$add_argument('-f', '--formula', default='~E1+E2+E3+E4+E5+E6',
                    help='Formula describing enhancer interactions. Default ~E1+E2+E3+E4+E5+E6.')
parser$add_argument('-ab', '--activitybounds',  default='c(-150, 150)',
                    help='Activity bounds to pass to superEnhancerDataObject() function. Default [-150, 150].')
parser$add_argument('-eb', '--errorbounds',  default='c(10^-3, 10^3)',
                    help='Error parameter bounds to pass to superEnhancerDataObject() function. Default [0.001, 1000].')
parser$add_argument('-sb', '--scalebounds',  default='c(10^-3, 10^3)',
                    help='Scale parameter bounds to pass to superEnhancerDataObject() function. Default [0.001, 1000].')
# get command line options, if help encountered print help and exit,
#   otherwise if options not found on command line, set defaults
args <- parser$parse_args()
# read data into df
datain <- read.csv(args$data)

r <- as.numeric(args$replicates)
wt <- args$wtnorm
i <- as.numeric(args$iterations)
f <- eval(parse(text=args$formula))
ab <- eval(parse(text=args$activitybounds))
eb <- eval(parse(text=args$errorbounds))
sb <- eval(parse(text=args$scalebounds))

# define parameters to test
error.models <- c('gaussian','lognormal')
link.functions <- c('additive','exponential','logistic')

# do some testing with contrived  Bglobin data:
test.Bglobin <- function(expr.reps, wt.norm, optim.iter, out = 'outputs/', enhancer.formula, activity.bounds = c(-150,150), error.bounds = c(10^-3, 10^3), scale.bounds = c(10^-3, 10^3)){
  # expr.reps = number of replicates of expression data for each enhancer condition. can be list. 
  # wt.norm = if TRUE, generate normally-distributed datapoints around 1 to represent WT expression. can be list. 
  # optim.iter = total number of iteration for optimization function to perform. can be list. 
  # out = path to output directory
  # enhancer.formula = interactions between variables to be modeled
  # maxit = maximum iterations of optimDE() to run
  # activity.bounds, error.bounds, scale.bounds = options to pass to enhancerDataObject() that restrict search for coefficients
  
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
        ggsave(figure, filename = paste0(out, 'BglobinSummary_', reps, 'reps_norm', norm.strategy, '_', iterations, 'iter_', Sys.Date(), '.pdf'), 
               device = 'pdf', width = 12, height = 8, units = 'in')
      }
    }
  }
  # append metadata to df for export
  master.out %>%
    mutate(enhancer.formula = deparse(enhancer.formula), activity.bounds = deparse(activity.bounds), error.bounds = deparse(error.bounds), scale.bounds = deparse(scale.bounds)) -> master.out
  # export data
  write.csv(master.out, file = paste0(out,'BglobinResults_',Sys.Date(),'.csv'), row.names = F)
}

test.Bglobin(expr.reps = r, wt.norm = wt, df = datain, optim.iter = i, out = args$outputlocation, 
             enhancer.formula = f, activity.bounds = ab, error.bounds = eb, scale.bounds = sb)
