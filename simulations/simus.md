---
title: "First test with graph.cpt package"
subtitle: "Case of change in a sequence of SPD matrices"
author: "Vincent Runge"
date: "09/05/2024"
output:
  html_document:
    keep_md: true
    css: styles.css
    toc: true
    toc_float: true
    highlight: tango
    number_sections: true
---





```r
library(graph.cpt) #our package
```

The outline of this document is as follows. In a first example, we test our functions on a simple case. 

- We generate data (time series) with a given SPD structure, with a few changes

- We define the graph structure used to generate the data (as if we knew the state values, oracle case)

- We use the `graph_cpt_manifold` function to infer the initial change point structure

- We discuss the result

In a second example, we also search a way to infer good state values for the graph. 

# A simple example

We test the `dust_R_1D` function with different pruning options and different data models.


```r
data <- dataGenerator_SPD(chpts = 10, d = 3, eigenvalues_mat = matrix(c(1,2,3),3,1))
data
```

```
## [[1]]
##           [,1]       [,2]       [,3]
## [1,] 2.5899534  0.1975633  0.5770345
## [2,] 0.1975633  1.2274095 -0.3726275
## [3,] 0.5770345 -0.3726275  2.1826371
## 
## [[2]]
##             [,1]        [,2]       [,3]
## [1,]  2.69456154 -0.05607395  0.4836073
## [2,] -0.05607395  1.28634102 -0.5166574
## [3,]  0.48360734 -0.51665738  2.0190974
## 
## [[3]]
##            [,1]       [,2]       [,3]
## [1,]  1.4929177 -0.1320074 -0.5031646
## [2,] -0.1320074  2.6214081 -0.6334224
## [3,] -0.5031646 -0.6334224  1.8856743
## 
## [[4]]
##            [,1]       [,2]       [,3]
## [1,]  2.3049207 -0.7289568 -0.3760932
## [2,] -0.7289568  2.2324585 -0.3305766
## [3,] -0.3760932 -0.3305766  1.4626208
## 
## [[5]]
##             [,1]       [,2]        [,3]
## [1,]  1.42990468  0.7993515 -0.01684437
## [2,]  0.79935146  2.5283032 -0.24046567
## [3,] -0.01684437 -0.2404657  2.04179211
## 
## [[6]]
##            [,1]      [,2]       [,3]
## [1,]  1.8004937 0.3505328 -0.6028182
## [2,]  0.3505328 1.5476518  0.4231655
## [3,] -0.6028182 0.4231655  2.6518545
## 
## [[7]]
##           [,1]       [,2]       [,3]
## [1,] 1.9897965  0.5231619  0.7957016
## [2,] 0.5231619  1.7234148 -0.1171966
## [3,] 0.7957016 -0.1171966  2.2867887
## 
## [[8]]
##             [,1]        [,2]       [,3]
## [1,]  1.82411383 -0.04991726  0.9369619
## [2,] -0.04991726  2.05024851 -0.3082008
## [3,]  0.93696194 -0.30820084  2.1256377
## 
## [[9]]
##           [,1]       [,2]       [,3]
## [1,] 2.7663596  0.3727490  0.4706976
## [2,] 0.3727490  2.0232345 -0.1839813
## [3,] 0.4706976 -0.1839813  1.2104059
## 
## [[10]]
##            [,1]       [,2]       [,3]
## [1,]  1.8386430  0.2210293 -0.9521002
## [2,]  0.2210293  2.0719799 -0.1582673
## [3,] -0.9521002 -0.1582673  2.0893770
```

