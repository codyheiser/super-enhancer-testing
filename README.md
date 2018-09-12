### super-enhancer-testing
Testing models of super enhancers using superEnhancerModelR package from Dukler, et al. (2017)  
  
##### `superE_datagen_CH.R`
Contains functions for generating expression data for testing the super enhancer model.  
Based on code from AP (2017), and beta-globin expression data from Bender, et al (2012).  
Also defines functions for reshaping data into format amenable to `superEnhancerModelR` package as well as personalized plotting functions adapted from the package itself.  

##### `superE_testing_CH.R`
Uses functions defined in `superE_datagen_CH.R` to test modeling of contrived beta-globin data.  
