### super-enhancer-testing
Testing models of super enhancer activity using superEnhancerModelR package from Dukler, et al. (2017)  
  
##### `utilityfunctions_superE.R`
* Contains functions for generating B-globin data for testing the super enhancer model based on code from AP (2017), and expression data from Bender, et al (2012)  
* Defines functions for reshaping data into format amenable to `superEnhancerModelR` package  
* Functions for streamlining model generation using `superEnhancerModelR`  
* Personalized plotting functions adapted from the package itself  
  
##### `test_superE.R`
Command-line generation of `superEnhancerModelR` models using all combinations of Link and Error functions, returning plots of model fits and BIC/relative BIC values and a _.csv_ file of coefficients and metadata for each model.  
  
usage:  
```
Rscript test_superE.R [-h] [-i ITERATIONS [ITERATIONS ...]] [-f FORMULA] [-ab ACTIVITYBOUNDS] [-eb ERRORBOUNDS] [-sb SCALEBOUNDS] data outputlocation
```
