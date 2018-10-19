# Compilation of results and figure generation for superEnhancerModelR package
# @author: C Heiser
# 07October2018

rm(list=ls())
source('utilityfunctions_superE.R')
require(tools)
require(plotly)

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

# split vector by delimeter, return element as vector
vectorsplit <- function(v, delim = '\\_', keep = 1){
  # v = vector of strings to split
  # delim = delimiter
  # keep = which element to return as list
  return(sapply(strsplit(as.character(v),delim), `[`, keep))
}

# split a column of strings from a data.frame to two columns by a delimiter
columnsplit <- function(df, clmn, newnames = c('name1', 'name2'), drop.orig = F, ...){
  # df = data.frame to operate on
  # clmn = name of column in df to split
  # newnames = vector containing new column names, in order
  # drop.orig = remove original clmn from data.frame?
  # ... = options to pass to vectorsplit function; specifically "delimiter = '\\_'"
  eval(parse(text = paste0('df$',newnames[1],'<-vectorsplit(df$',clmn,',...,keep=1)'))) # put string from before delim into column with name1
  eval(parse(text = paste0('df$',newnames[2],'<-vectorsplit(df$',clmn,',...,keep=2)'))) # put string from after delim into column with name2
  if(drop.orig){
    eval(parse(text = paste0('df<-subset(df, select = -',clmn,')'))) # drop original column from df if drop.orig flag set to TRUE
  }
  return(df)
}

# plot raw expression values from df in superE format 
plot.raw.expr <- function(superEdata, plot.title){
  superEdata %>%
    group_by(condition) %>%
    summarise(expr.mean = mean(expression), expr.se = sd(expression)/sqrt(n()), n = n()) -> superEsum
  return(
    ggplot(data = superEsum, aes(x = condition, y = expr.mean))+
      geom_jitter(data = superEdata, aes(x = condition, y = expression), width = 0.15, size = 2.5, alpha = 0.6, color = 'goldenrod')+
      geom_bar(color='black', alpha = 0, width = 0.6, stat = 'identity')+
      geom_errorbar(color='black', width = 0.25, aes(ymin = expr.mean-expr.se, ymax = expr.mean+expr.se))+
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
read.all(filetype = 'csv',
         dir = '~/Dropbox/_Venters_Lab_Resources/3_Rotation_Students/4_Cody/superE/Bglobin_30Sep18/') %>%
  rename(`x-Int.`=X.Intercept., Scale=scale, `Activity Bounds`=activity.bounds) %>% # clean up some column names
  mutate(temp = vectorsplit(vectorsplit(file, "\\."), "X", 2)) %>% # extract replicate and WT normalization from filename
  columnsplit(clmn = 'temp', newnames = c('Replicates','WT'), drop.orig = T, delim = "\\_") %>%
  mutate(act = as.numeric(vectorsplit(vectorsplit(`Activity Bounds`, ',', 2), '\\)')))-> master

# determine lowest BIC in each test cohort and create 'winner' column in master df
master %>% 
  group_by(file, optim.iter) %>% 
  filter(bic==min(bic)) %>%
  mutate(winner = T) -> low.bic
master %>%
  group_by(file, optim.iter) %>%
  filter(bic!=min(bic)) %>%
  mutate(winner = F) %>%
  bind_rows(low.bic) -> master

#########################################################################################################################
# try some plotting stuff
p <- plot_ly(master, x = ~optim.iter, y = ~act, z = ~bic, color = ~winner) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'DEoptim Iterations'),
                      yaxis = list(title = 'Activity Bounds'),
                      zaxis = list(title = 'BIC')))

master %>% 
  group_by(act, optim.iter) %>% 
  filter(bic==max(bic)) %>%
  select(bic, optim.iter) %>%
  spread(value = bic, key = optim.iter) %>%
  ungroup() %>%
  select(-act) -> surface.hi

master %>% 
  group_by(act, optim.iter) %>% 
  filter(bic==min(bic)) %>%
  select(bic, optim.iter) %>%
  spread(value = bic, key = optim.iter) %>% 
  ungroup() %>%
  select(-act) -> surface.lo

# 3D surface of best and worst BICs for each condition across all tests
p <- plot_ly(x = unique(master$act), y = unique(master$optim.iter)) %>% 
  add_surface(z = as.matrix(surface.lo), showscale=F) %>%
  add_surface(z = as.matrix(surface.hi), showscale=F) %>%
  layout(showlegend=F, scene = list(xaxis = list(title = 'Activity Bounds'),
                      yaxis = list(title = 'Optim. Iterations'),
                      zaxis = list(title = 'BIC')))

# heatmap?
row.names(surface.hi) <- unique(master$act)
heatmap(as.matrix(surface.hi), xlab="Iterations", ylab="Bounds", Rowv = NA, Colv = NA)

# heatmap?
row.names(surface.lo) <- unique(master$act)
heatmap(as.matrix(surface.lo), xlab="Iterations", ylab="Bounds", Rowv = NA, Colv = NA)

#########################################################################################################################
# write false positive BIC calls to .csv file
low.bic %>%
  select(link, error, optim.iter, bic, rel.bic, `Activity Bounds`, Replicates, WT) %>%
  filter(link!='additive') %>%
  write.csv('outputs/false_positives.csv', row.names = F)

# generate plot of BIC values
master %>%
  filter(error=='lognormal') %>% # ignore gaussian error for plotting simplicity
  ggplot(aes(x = link, y = bic, color = `Activity Bounds`, shape = winner))+
  scale_color_manual(values = c('c(-10, 10)'='firebrick1', 'c(-50, 50)'='firebrick3', 'c(-100, 100)'='goldenrod1',
                     'c(-150, 150)'='goldenrod2', 'c(-500, 500)'='goldenrod3', 'c(-1000, 1000)'='goldenrod4'))+
  geom_jitter(width = 0.2, size = 3, alpha = 0.6)+
  labs(x = NULL, y = 'BIC', shape = 'Best Fit', title = 'BIC Values for Log-Normal Error')+
  plot.opts+
  theme(axis.text.x = element_text(angle = 50, hjust = 1), 
        panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank())+
  facet_grid(Replicates~optim.iter, scales = 'free') -> lognormal.bic.plt
ggsave(plot = lognormal.bic.plt, filename = 'outputs/lognormal_bic.pdf', device = 'pdf', height = 6, width = 8, units = 'in')

# generate summary figure of link coefficients to show convergence based on replicates, bounds, and iterations
plt.list <- list()
for(bounds in unique(master$`Activity Bounds`)){
  master %>%
    filter(`Activity Bounds`==bounds) %>%
    gather(key = 'variable', value = 'value', `x-Int.`, E1, E2, E3, E4, E5, E6, Scale) %>%
    ggplot(aes(x = variable, y = value, color = optim.iter))+
    geom_jitter(width = 0.2,size = 2.5, alpha = 0.6)+
    facet_grid(link~Replicates, scales = 'free')+
    labs(x=NULL,y=NULL, title = bounds, color = 'Optimization\nIterations')+
    plot.opts +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank()) -> plt
  assign(paste0(bounds,'.plt'), plt)
  plt.list <- append(plt.list, paste0(bounds,'.plt'))
}
coeff.fig <- ggarrange(plotlist = eval(parse(text = paste0('list(`',paste(unlist(plt.list), collapse = '`,`'),'`)'))), 
                     ncol = 2, nrow = 3, align = "v", labels = 'auto', common.legend = T, legend = 'right') %>%
  annotate_figure(left = text_grob('Value', rot = 90, size = 14),
                  bottom = text_grob('Link Coefficient', size = 14))
ggsave(plot = coeff.fig, filename = 'outputs/link_coeff_fig.pdf', device = 'pdf', height = 14, width = 16, units = 'in')

#########################################################################################################################
# reshape some novel enhancer data for testing
slc25a37 <- read.default('inputs/Slc25a37-enhDeletion_modCH.csv')

slc25a37 %>%
  select(-grep(pattern = 'mean|sd|sem', names(slc25a37))) %>%
  gather(key = variable, value = expression, -TimePoint) %>%
  filter(!is.na(expression)) %>%
  mutate_cond(str_detect(variable, 'wt'), variable = 'wt') %>%
  mutate(condition = vectorsplit(variable, '\\.')) %>%
  mutate(E1=1, E2=1, E3=1, P=1) %>%
  mutate_cond(str_detect(condition,'1'), E1=0) %>%
  mutate_cond(str_detect(condition,'2'), E2=0) %>%
  mutate_cond(str_detect(condition,'3'), E3=0) %>%
  mutate_cond(str_detect(condition,'promoter'), P=0) %>%
  select(TimePoint, expression, condition, E1, E2, E3, P) -> slc25a37_superE

slc25a37_superE %>%
  filter(TimePoint == '0h') %>%
  select(-TimePoint) %>%
  write.csv('inputs/slc25a37_0hr.csv', row.names = F)
slc25a37_superE %>%
  filter(TimePoint == '4h') %>%
  select(-TimePoint) %>%
  write.csv('inputs/slc25a37_4hr.csv', row.names = F)
slc25a37_superE %>%
  filter(TimePoint == '8h') %>%
  select(-TimePoint) %>%
  write.csv('inputs/slc25a37_8hr.csv', row.names = F)
slc25a37_superE %>%
  filter(TimePoint == '12h') %>%
  select(-TimePoint) %>%
  write.csv('inputs/slc25a37_12hr.csv', row.names = F)
slc25a37_superE %>%
  filter(TimePoint == '24h') %>%
  select(-TimePoint) %>%
  write.csv('inputs/slc25a37_24hr.csv', row.names = F)
slc25a37_superE %>%
  filter(TimePoint == '48h') %>%
  select(-TimePoint) %>%
  write.csv('inputs/slc25a37_48hr.csv', row.names = F)

slc25a37_superE$TimePoint = factor(slc25a37_superE$TimePoint, levels=c('0h','4h','8h','12h','24h','48h'))
slc25a37_superE %>%
  group_by(condition, TimePoint) %>%
  summarise(expr.mean = mean(expression), expr.se = sd(expression)/sqrt(n()), n = n()) %>%
  ggplot(aes(x = condition, y = expr.mean))+
  geom_jitter(data = slc25a37_superE, aes(x = condition, y = expression), width = 0.15, size = 2.5, alpha = 0.6, color = 'goldenrod')+
  geom_bar(color='black', alpha = 0, width = 0.6, stat = 'identity')+
  geom_errorbar(color='black', width = 0.25, aes(ymin = expr.mean-expr.se, ymax = expr.mean+expr.se))+
  labs(x = NULL, y = 'Expression', title = expression(italic('slc25a37')))+
  facet_wrap(~TimePoint, ncol = 3, scales = 'free_y')+
  plot.opts+
  theme(axis.text.x = element_text(angle=50, hjust=1)) -> slc25a37.plt.timepoint
ggsave(plot = slc25a37.plt.timepoint, filename = 'outputs/slc25a37_timepointplot_raw.pdf', device = 'pdf', width = 8, height = 6, units = 'in')

slc25a37_superE %>%
  group_by(condition) %>%
  summarise(expr.mean = mean(expression), expr.se = sd(expression)/sqrt(n()), n = n()) -> superEsum
ggplot(data = superEsum, aes(x = condition, y = expr.mean, color = TimePoint))+
  geom_jitter(data = slc25a37_superE, aes(x = condition, y = expression), width = 0.15, size = 2.5, alpha = 0.6)+
  geom_bar(color='black', alpha = 0, width = 0.6, stat = 'identity')+
  geom_errorbar(color='black', width = 0.25, aes(ymin = expr.mean-expr.se, ymax = expr.mean+expr.se))+
  labs(x = NULL, y = 'Expression', title = expression(italic('slc25a37')))+
  plot.opts+
  theme(axis.text.x = element_text(angle=50, hjust=1)) -> slc25a37.plt.all
ggsave(plot = slc25a37.plt.all, filename = 'outputs/slc25a37_raw.pdf', device = 'pdf', width = 6, height = 4, units = 'in')
