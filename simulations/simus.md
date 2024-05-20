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

We generate 10 matrices of size 3, with eigenvalues 1,2,3.


```r
data <- dataGenerator_SPD(chpts = 10, d = 3, eigenvalues_mat = matrix(c(1,2,3),3,1))
data
```

```
## [[1]]
##           [,1]       [,2]       [,3]
## [1,] 1.8132507  0.2477274  0.3008909
## [2,] 0.2477274  2.2869722 -0.8856971
## [3,] 0.3008909 -0.8856971  1.8997771
## 
## [[2]]
##            [,1]       [,2]        [,3]
## [1,] 1.52616289  0.8735854  0.02380141
## [2,] 0.87358538  2.4729052 -0.11045889
## [3,] 0.02380141 -0.1104589  2.00093191
## 
## [[3]]
##            [,1]      [,2]       [,3]
## [1,]  2.1513670 0.1917570 -0.8849441
## [2,]  0.1917570 1.8008347  0.3842805
## [3,] -0.8849441 0.3842805  2.0477983
## 
## [[4]]
##            [,1]       [,2]       [,3]
## [1,]  1.9224214  0.1772623 -0.6789586
## [2,]  0.1772623  2.4200540 -0.5980942
## [3,] -0.6789586 -0.5980942  1.6575246
## 
## [[5]]
##            [,1]      [,2]       [,3]
## [1,]  2.0487732 0.4916912 -0.3181710
## [2,]  0.4916912 2.7394974  0.2677761
## [3,] -0.3181710 0.2677761  1.2117294
## 
## [[6]]
##            [,1]       [,2]       [,3]
## [1,]  2.1228155 -0.4512562 -0.1225418
## [2,] -0.4512562  1.8542040  0.8734514
## [3,] -0.1225418  0.8734514  2.0229805
## 
## [[7]]
##           [,1]       [,2]       [,3]
## [1,] 2.0678286  0.4428757  0.8168337
## [2,] 0.4428757  1.6737256 -0.2184462
## [3,] 0.8168337 -0.2184462  2.2584458
## 
## [[8]]
##            [,1]       [,2]       [,3]
## [1,]  1.9943567 -0.5172189 -0.7576916
## [2,] -0.5172189  1.6315774 -0.1433243
## [3,] -0.7576916 -0.1433243  2.3740660
## 
## [[9]]
##           [,1]      [,2]      [,3]
## [1,] 1.7575524 0.7269801 0.2087124
## [2,] 0.7269801 1.7632347 0.5056651
## [3,] 0.2087124 0.5056651 2.4792129
## 
## [[10]]
##            [,1]       [,2]      [,3]
## [1,]  1.7247008 -0.4758645 0.4229626
## [2,] -0.4758645  2.8088267 0.2955243
## [3,]  0.4229626  0.2955243 1.4664725
```
We compute the distances to reference matrices `sigma1`, `sigma2`, `sigma3`.


```r
dim <- 3
s1 <- 1*sample(dim)
s2 <- 1.5*sample(dim)
s3 <- 2*sample(dim)

s <- matrix(c(s1,s2,s3),dim,3)

sigma1 <- diag(s1, dim)
sigma2 <- diag(s2, dim)
sigma3 <- diag(s3, dim)
```


The distances are thus:


```r
#for(i in 1:10){print(eigen(data[[i]]))}
sigma123 <- list(sigma1,sigma2,sigma3)
dists <- SPDts_to_dists(data, sigma123)
dists
```

```
##          [,1]     [,2]     [,3]      [,4]      [,5]     [,6]     [,7]     [,8]
## [1,] 1.032169 1.095759 1.133272 0.8673264 0.5599761 1.174373 1.240442 1.288210
## [2,] 1.417830 1.401100 1.301320 1.5078121 1.6214078 1.267316 1.207422 1.163191
## [3,] 1.595983 1.457891 1.678272 1.5705534 1.5764515 1.724805 1.681154 1.674893
##          [,9]     [,10]
## [1,] 1.329925 0.7358651
## [2,] 1.144344 1.5876274
## [3,] 1.631696 1.4825168
```

```r
apply(dists, 2, which.min)
```

```
##  [1] 1 1 1 1 1 1 2 2 2 1
```





