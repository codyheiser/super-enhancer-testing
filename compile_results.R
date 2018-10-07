# Compilation of results and figure generation for superEnhancerModelR package
# @author: C Heiser
# 07October2018

rm(list=ls())
source('utilityfunctions_superE.R')

#########################################################################################################################
# Define functions

# read in .csv or .xlsx file
read.default <- function(file, ...){
  # file = Path and filename plus extension. May be .csv or .xlsx
  # sheet = for Excel files, name or index of sheet to read
  
  # timer
  ptm <- proc.time()
  
  if(file_ext(file) == 'csv'){
    df <- read.csv(file, check.names = T, stringsAsFactors = F, ...)
  }else{if(file_ext(file) == 'txt'){
    df <- read.table(file, header = T, sep = '\t', ...)
  }else{if(file_ext(file) %in% c('xls', 'xlsx')){
    df <- read_excel(file, ...)
  }else{
    df <- NA
  }}}
  print(proc.time() - ptm) # see how long this took
  return(df)
}

# read in all files of common type from folder, concatenate by row
read.all <- function(filetype, dir = '.', ...){
  # filetype = name of extension to read ('csv', 'xls', or 'xlsx'). you can also use globs (e.g. 'myfile*.csv')
  # dir = directory to read files from
  
  # timer
  ptm <- proc.time()
  
  vars <- list() # initiate empty list for later concatenation of dfs
  for(f in list.files(path = dir, pattern = filetype, recursive = T)){
    if(str_detect(f, regex('\\/\\~\\$'))){
      # ignore files that are open by Windows
    }else{
      name <- make.names(f) # get syntactically valid name of file
      print(paste0('Reading ',name))
      df <- read.default(file_path_as_absolute(paste0(dir,f)), ...) # read csv or Excel file into dataframe
      df$file <- name # create 'file' column that has metadata pointing to file name
      assign(name, df) # rename the df as the file ID
      vars <- append(vars, name) # add name of new df to list of variables
    }
  }
  # concatenate rows of all dfs
  combined <- eval(parse(text = paste0('do.call(rbind.fill, args = list(',paste(vars,collapse = ','),'))')))
  print(proc.time() - ptm) # see how long this took
  return(combined)
}

# plot raw expression values from df in superE format 
plot.raw.expr <- function(superEdata, plot.title){
  return(
    ggplot(data = superEdata, aes(x = condition, y = expression))+
      geom_jitter(width = 0.15, size = 2.5, alpha = 0.6, color = 'goldenrod')+
      geom_boxplot(color='black', alpha = 0, width = 0.5, outlier.stroke = 0, outlier.size = 0)+
      labs(x = NULL, y = 'Expression', title = plot.title)+
      plot.opts
  )
}

#########################################################################################################################
# look at raw expression data
# read in data
wap <- read.csv('inputs/wap_data_superEformat.csv')
a.globin <- read.csv('inputs/alpha_globin_data_superEformat.csv')
b.globin.10 <- read.csv('inputs/beta_globin_10_norm.csv')
b.globin.100 <- read.csv('inputs/beta_globin_100_norm.csv')
# generate plots
wap.plt <- plot.raw.expr(wap, expression(italic('wap')))
a.globin.plt <- plot.raw.expr(a.globin, expression(paste(alpha,'-globin')))
b.globin.plt.1 <- plot.raw.expr(b.globin.10, expression(paste(beta,'-globin (n=10)')))
b.globin.plt.2 <- plot.raw.expr(b.globin.100, expression(paste(beta,'-globin (n=100)')))
# compile figure
raw.fig <- ggarrange(plotlist = lapply(list(wap.plt, a.globin.plt, b.globin.plt.1, b.globin.plt.2), 
                                       FUN = function(x){return(x+labs(y=NULL)+theme_pubr())}), 
                     ncol = 2, nrow = 2, align = "v", labels = c('A','B','C','D')) %>%
  annotate_figure(left = text_grob('Expression', rot = 90, size = 14),
                  bottom = text_grob('Enhancer Status', size = 14))
ggsave(plot = raw.fig, filename = 'outputs/raw_expression_fig.pdf', device = 'pdf', height = 6, width = 8, units = 'in')

#########################################################################################################################
# compile results from b.globin testing 
master <- read.all(filetype = 'csv', 
                   dir = '~/Dropbox/_Venters_Lab_Resources/3_Rotation_Students/4_Cody/superE/Bglobin_30Sep18/')

plt.list <- list()
for(bounds in unique(master$activity.bounds)){
  master %>%
    rename(`x-Int.`=X.Intercept., Scale=scale, `Activity Bounds`=activity.bounds) %>%
    filter(`Activity Bounds`==bounds) %>%
    gather(key = 'variable', value = 'value', `x-Int.`, E1, E2, E3, E4, E5, E6, Scale) %>%
    ggplot(aes(x = variable, y = value, color = optim.iter))+
    geom_jitter(width = 0.2,size = 2.5, alpha = 0.6)+
    facet_grid(link~error, scales = 'free')+
    labs(x=NULL,y=NULL, title = bounds)+
    plot.opts +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) -> plt
  assign(paste0(bounds,'.plt'), plt)
  plt.list <- append(plt.list, paste0(bounds,'.plt'))
}
coeff.fig <- ggarrange(plotlist = eval(parse(text = paste0('list(`',paste(unlist(plt.list), collapse = '`,`'),'`)'))), 
                     ncol = 2, nrow = 3, align = "v", labels = 'auto', common.legend = T, legend = 'right') %>%
  annotate_figure(left = text_grob('Value', rot = 90, size = 14),
                  bottom = text_grob('Link Coefficient', size = 14))
ggsave(plot = coeff.fig, filename = 'outputs/link_coeff_fig.pdf', device = 'pdf', height = 16, width = 18, units = 'in')



