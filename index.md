---
title       : 140.655 Biostat Lab 5
author      : L. Collado-Torres
framework   : minimal
highlighter : prettify
hitheme     : twitter-bootstrap
mode        : selfcontained
github      : {user: lcolladotor, repo: BiostatLab5, branch: gh-pages}
widgets     : [mathjax, disqus, ganalytics]
assets:
  css: 
    - "http://fonts.googleapis.com/css?family=PT+Sans"
    - "http://odyniec.net/articles/turning-lists-into-trees/css/tree.css"
---













This project was done by [Emily Huang](http://www.jhsph.edu/departments/biostatistics/directory/students/phd.html) and [Leonardo Collado-Torres](http://bit.ly/LColladoTorres) for the [140.655 Analysis of Longitudinal Data](http://www.jhsph.edu/courses/course/140.655/01/2013/17988/) laboratory session 5 for biostatistics students.


## Instructions

```
Conduct, using a simulation, a sensitivity analysis of the change in effects
estimates as a function of the number of quadrature points.
Consider different integration methods.
```

As a reference/dataset starting point consider using (<span class="showtooltip" title="Lesaffre E and Spiessens B (2001). 'on The Effect of The Number of Quadrature Points in A Logistic Random Effects Model: an Example.' Journal of The Royal Statistical Society: Series C (Applied Statistics), 50, pp. 325-335. ISSN 0035-9254."><a href="http://dx.doi.org/10.1111/1467-9876.00237">Lesaffre & Spiessens, 2001</a></span>).


## Theory

Below some of the key theoretical points are shown. For more details, please check out the [presentation](http://lcolladotor.github.io/BiostatLab5/GLMM.pdf).

### General GLMM setup

1) Distribution assumption
* $Y_{ij}|b_{i}\sim$ Exponential family
* $Var(Y_{ij}|b_{i})=v\{E(Y_{ij}|b_{i})\}\phi$, where $v$ is a known function
* $Cov(Y_{ij},Y_{ik}|b_{i})=0$

2) Systematic component
* $\eta_{ij}=X_{ij}\beta+Z_{ij}b_{i}$

3) Link function
* $g\{E(Y_{ij}|b_{i})\}=\eta_{ij}=X_{ij}\beta+Z_{ij}b_{i}$ for some known link function, $g$

4) Random effects
* Assumed to have some probability distribution, such as $b_{i}\sim MVN(0,G)$
* $b_{i}$ are assumed to be independent of the covariates


### Gauss-Hermite Quadrature

$$\int_{-\infty}^{\infty} & h(v)e^{-v^{2}}dv & \approx \sum_{k=1}^{d}h(x_{k})w_{k}$$

* $d$ quadrature points (weights, $w_{k}$, and evaluation points, $x_{k})$
* The more quadrature points used, the more accurate the approximation
* But computational burden increases with quadrature points, and grows exponentially with the number of random effects


## Simulation



## References

Web document generated using `slidify` (<span class="showtooltip" title="Vaidyanathan R (2012). slidify: Generate reproducible html5 slides from R markdown. R package version 0.4."><a href="http://ramnathv.github.com/slidify/">Vaidyanathan, 2012</a></span>). Citations made with `knitcitations` (<span class="showtooltip" title="Boettiger C (2014). knitcitations: Citations for knitr markdown files. R package version 0.5-0."><a href="http://CRAN.R-project.org/package=knitcitations">Boettiger, 2014</a></span>). 



- Carl Boettiger,   (2014) knitcitations: Citations for knitr markdown files.  [http://CRAN.R-project.org/package=knitcitations](http://CRAN.R-project.org/package=knitcitations)
- Emmanuel Lesaffre, Bart Spiessens,   (2001) on The Effect of The Number of Quadrature Points in A Logistic Random Effects Model: an Example.  *Journal of The Royal Statistical Society: Series C (Applied Statistics)*  **50**  325-335  [10.1111/1467-9876.00237](http://dx.doi.org/10.1111/1467-9876.00237)
- Ramnath Vaidyanathan,   (2012) slidify: Generate reproducible html5 slides from R markdown.  [http://ramnathv.github.com/slidify/](http://ramnathv.github.com/slidify/)




## R code




## Reproducibility

This report was last updated on


```
## [1] "2014-03-09 20:46:21 EDT"
```


R session information:


```
## R version 3.0.2 (2013-09-25)
## Platform: x86_64-apple-darwin10.8.0 (64-bit)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] knitcitations_0.5-0 bibtex_0.3-6        knitr_1.5          
## [4] slidify_0.4         lme4_1.0-6          Matrix_1.1-2-2     
## [7] lattice_0.20-23    
## 
## loaded via a namespace (and not attached):
##  [1] codetools_0.2-8 digest_0.6.4    evaluate_0.5.1  formatR_0.10   
##  [5] grid_3.0.2      httr_0.2        markdown_0.6.4  MASS_7.3-29    
##  [9] minqa_1.2.3     nlme_3.1-111    Rcpp_0.11.0     RCurl_1.95-4.1 
## [13] splines_3.0.2   stringr_0.6.2   tools_3.0.2     whisker_0.3-2  
## [17] XML_3.95-0.2    xtable_1.7-1    yaml_2.1.10
```



```r
## Generate this report
library("slidify")
slidify("index.Rmd")

```



<div id='disqus_thread'></div>


