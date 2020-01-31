# CircMLE

Maximum Likelihood Analysis of Circular Data

## Description
A series of wrapper functions to implement the 10 maximum likelihood models of animal orientation described by Schnute and Groot (1992) doi: [10.1016/S0003-3472(05)80068-5](https://doi.org/10.1016/S0003-3472(05)80068-5). The functions also include the ability to use different optimizer methods and calculate various model selection metrics (i.e., AIC, AICc, BIC).
This framework is designed for modeling any dataset represented by angles (e.g, orientation, periodic, etc) using the above models. Main features are listed as follows.

- Calculate the likelihood of any one or all of the 10 models of orientation
- Compare any two nested models using a likelihood ratio test
- Plot the observed dataset and any of the model-fitted results
- Calculate the Hermans-Rasson test or Pycke test for directionality


## Install CircMLE (from an R console)
- To install from CRAN
  * First install the R package 'circular' from CRAN using the command `install.packages("circular")`
  * Then install the CircMLE package using `install.packages("CircMLE")`
  * Load the package into your working R environment using `library(CircMLE)`


## Version History
- Version 0.2.3 2020/1/29
  * Added the ability to perform the Hermans-Rasson and Pycke tests using code kindly provided by Lukas Landler, Graeme Ruxton, and E. Pascal Malkemper.

- Version 0.2.2 2019/10/17
  * Improved communication between *CircMLE* and R 'circular' objects, especially for improved plotting when using 'template = "geographics"'.
  
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

<b><i>If using the Hermans-Rasson or Pycke tests then cite:</b></i>  
Landler, L., Ruxton, G. D., and Malkemper, E. P. (2019) The Hermans–Rasson test as a powerful alternative to the Rayleigh test for circular statistics in biology. BMC Ecology 19: 30; doi: [10.1186/s12898-019-0246-8](https://doi.org/10.1186/s12898-019-0246-8)

- Or enter the command `citation("CircMLE")` into your R console

## Contact
Robert Fitak  
Department of Biology  
University of Central Florida  
USA  
rfitak9@gmail.com  