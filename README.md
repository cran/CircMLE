# CircMLE

Maximum Likelihood Analysis of Circular Data

## Description
A series of wrapper functions to implement the 10 maximum likelihood models of animal orientation described by Schnute and Groot (1992) doi: [10.1016/S0003-3472(05)80068-5](https://doi.org/10.1016/S0003-3472(05)80068-5). The functions also include the ability to use different optimizer methods and calculate various model selection metrics (i.e., AIC, AICc, BIC).
This framework is designed for modeling any dataset represented by angles (e.g, orientation, periodic, etc) using the above models. Main features are listed as follows.

- Calculate the likelihood of any one or all of the 10 models of orientation
- Compare any two nested models using a likelihood ratio test
- Plot the observed dataset and any of the model-fitted results


## Install CircMLE (from an R console)
- To install from CRAN
  * First install the R package 'circular' from CRAN using the command `install.packages("circular")`
  * Then install the CircMLE package using `install.packages("CircMLE")`
  * Load the package into your working R environment using `library(CircMLE)`


## Version History
- Version 0.2.1 2018/02/20
  * Added support for data vectors with the "geographics" template set when plotting the modeled results.
  * Added publication information
  * Added the README.md file

- Version 0.2.0 2017/06/29
  * Added a plotting function to visualize the observed and modeled results

- Version 0.1, 2017/05/13
  * Released the first version


## Citation
Fitak, R. R. and Johnsen, S. (2017) Bringing the analysis of animal orientation data full circle: model-based approaches with maximum likelihood. Journal of Experimental Biology 220: 3878-3882; doi: [10.1242/jeb.167056](https://doi.org/10.1242/jeb.167056)

- Or enter the command `citation("CircMLE")` into your R console

## Contact
Robert Fitak  
Department of Biology  
Duke University  
USA  
rfitak9@gmail.com  