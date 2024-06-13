###### TESTING THE CODE FROM RYAN, KILLICK - how are we doing it ######

library(devtools)
install_github("s-ryan1/Covariance_RMT_simulations/cpt.cov")

#THINGS TO CONSIDER:

#### 1. 
#They have two different settings simulating the matrix sequences
#one to satisfy the conditions from Wang et al, the second one to satisfy the conditions by
#Ryan, Killick


#### 2.
#The data used in the paper is n*p so we have n p-dimensional observations
#and the covariance matrices are calculated then from the original observations as part of the
#binary segmentation algorithm. When an interval (s,e) is tested, the covariance matrices are computed
#on (s,b), (b,e) and their statistic is a function of the two covariance matrices estimates.

#So not using a sliding window like we do (for applications)
#After generating the sequence of covariance matrices in this way
#they do not actually correspond to our model - it's not a piecewise constant jump in covariance matrices


#### 3.
#See these files on GitHub:
#wang_multi_cpts.R
#ratio_multi_cpts.R
#These functions only deal with 4 CPs, I will just modify the part of their run.sim function 
#that generates the data, so that the pattern is repeated and we don't need to have just 4 CPs in a pattern

#### 4.
#Minimum segment length, used in function gen.cpts.locs, set to 30
#They consider p (matrix dim, 3 in our case) to be increasing with n for their theory
#they say that max(4p, 30) is the default value for the minseglen
#If minseglen is smaller, the FPR is very large.
#
#Some questions:
#How long should the intervals between the change-points be? 
#What width of the sliding window should we use?


#### 5.
#The criteria for comparing is
#TPR - number of correctly detected / total number of cps
#FPR - number of falsely detected / total number of detected cps
#MAE - mean absolute error - the mean distance between all the estimated vs true covariance matrices

########## QUICK TEST - COPIED ratio_example.R FILE ########

library(cpt.cov)

bonferoni = function(n,alpha){ qnorm(1-.05/n)}

n=500
p=30
Y = matrix(rnorm(n*p,0,1), nrow=n)
Y[300:500,] = 2*Y[300:500,]

minseglen = 120
result =  matrix.dist.test.stat(Y, minseglen-p)

if(max(result, na.rm=T)>bonferoni(n,.05)){
  print(sprintf("Change detected at %d", which.max(abs(result))))
} else {
  print(sprintf("No change detected"))
} 

##### MODIFYING ratio_multi_cps.R  #####

library(cpt.cov)
library(tidyverse)

#args = commandArgs(trailingOnly=TRUE)
#iter = as.numeric(args[1])

#Compared to their run.sim function I added the arguments
#1. cp_num - number of states - 1
#2. num_repeats - number of times we are repeating the pattern

#this function generates the matrices that are far enough apart with respect to d2 
#(the metric proposed by ryan killick)

run.sim = function(i, n, p, cp_num, num_repeats, diff){
  
  if(n < 5*p*log(n)){return(NA) }
  sigma=list()
  
 
  # B is the matrix of eigenvectors of a pxp covariance matrix (random orthonormal matrix) 
  B = rnorm(10*p,0,1) %>% matrix(ncol=p) %>% cov %>% eigen %>% pluck("vectors")	
  #diagonal eigenvalue matrix (eigenvalues from U[0.1,10])
  sigma[[1]] = diag(runif(p,.1,10))  
  #sigma[[1]] is Lambda_1 in the text (in the code they later multiply with B to get Sigma_1 from the text)
  
  #The next sigma (Lambda) depends on the previous to control that the size of the changes is big enough
  #The changes are big enough if the eigenvalues are sufficiently far apart
  #diff = 1 #this is K (kappa) in the text?

  for(k in 2:(cp_num + 1)){
    delta = runif(p, 0, diff*p) %>% sort %>% diff
    delta = c(delta, diff*p-sum(delta)) + 1 #why do we have +1 here? Can't see it in the text
    delta = delta^((-1)^(rbernoulli(p,.5)))
    sigma[[k]] = delta*sigma[[k-1]]
  }
  #sigma up to here is the sequence of Lambda's
  sigma = map(sigma, ~B%*%((.x))%*%t(B)) #this is step 5 (page 32)
  sigma = rep(sigma, num_repeats)
  
  #gen.cpt.locs calls segment.remover which is in data_generation.R
  #segment.remover randomly picks a position i and removes an interval 
  #i-minseglen, i+minseglen from the sequence
  #this function is supposed to generate the change-point locations
  #but in a bit weird way?
  
  cpts =  gen.cpt.locs(n, cp_num, max(30,floor(p*log(n))))
  cpts = c(cpts, n)
  #i have to modify this, cpts contains the change-points of the first period only
  #have to add the change-points in other periods
  cpts_complete = cpts
  if(num_repeats > 1)
    for(i in 1:(num_repeats-1))
      cpts_complete = append(cpts_complete, n*i+cpts)
  cpts = cpts_complete[-length(cpts_complete)]
  
  lengths= diff(c(0,cpts,n*num_repeats))  #lengths of segments between cps
  
  #for each segment of length s, generate s*p normal variables, put them into matrix
  #each row is a N(0, I) random vector that we multiply by the corresponding Sigma to get the data as N(0, Sigma)
  noise = map(lengths, ~rnorm(.x*p, 0, 1)) %>% map(matrix, ncol=p)
  data = map2(noise, sigma, ~(.x)%*%.y) %>% reduce(rbind) 

  n = n*num_repeats  
  
  result.fisher = bin.seg(data, c(0,n), matrix.dist.test.stat, threshold=qnorm(1-.05/(n^2)), minseglen=max(c(30,4*p)), c())
  fisher = result.fisher %>% bin.seg.to.cpt(qnorm(1-.05/(n^2)))
  wang.thresh = wang.threshold(data)
  result.wang = bin.seg(data, c(0,n), wang.stat, threshold=wang.thresh, minseglen=p*log(n), c())
  wang = result.wang %>% bin.seg.to.cpt(wang.thresh)
  #print(fisher$cpts)
  result.galeano = bin.seg(data, c(0,n), cpt.cov:::galeano.cusum.stat, threshold=qnorm(1-.05/(n^2)), minseglen=40, c())
  galeano = result.galeano %>% bin.seg.to.cpt(qnorm(1-.05/(n^2)))
  
  
  fisher.cpt.error = detection.rates(fisher$cpts, cpts, 20)
  wang.cpt.error = detection.rates(wang$cpts, cpts, 20)
  galeano.cpt.error = detection.rates(galeano$cpts, cpts, 20)
  #set.seed(353*n + 541*p + 7*i)
  #return(max(matrix.dist.test.stat(data, 4*p), na.rm=T))
  
  if(p<=20 & p*(p+1) < 2*n){
    safe.bin.seg = safely(bin.seg)
    result.aue = safe.bin.seg(data, c(0,n), aue.stat, threshold=qnorm(.95), minseglen=p*log(n), c())
    if(length(result.aue[[1]])>0){
      aue = result.aue[[1]] %>% bin.seg.to.cpt(qnorm(.95))
      aue.cpt.error = detection.rates(aue$cpts, cpts, 20)
    }
    models = list(fisher$cpts, wang$cpts, aue$cpts, galeano$cpts)
    rates = map(models, detection.rates, cpts, 20)
    m_hat = map_dbl(rates, ~.x$m)
    TDR = map_dbl(rates, ~.x$TDR)
    FDR = map_dbl(rates, ~.x$FDR)
    mae = map_dbl(models, MAE,  cpts, data) 
    smae = map_dbl(models, spectral.error,  cpts, data) 
    result = tibble(
      method=c("Ratio", "Wang", "Aue", "Galeano"),
      TDR=TDR, FDR=FDR, m_hat=m_hat, 
      MAE=mae, SMAE=smae, n=n, p=p, power=diff)
    #result = result %>% gather("metric", "value",-method,-n,-p,-power)
    return(result)
  }
  
  models = list(fisher$cpts, wang$cpts, galeano$cpts)
  rates = map(models, detection.rates, cpts, 20)
  m_hat = map_dbl(rates, ~.x$m)
  TDR = map_dbl(rates, ~.x$TDR)
  FDR = map_dbl(rates, ~.x$FDR)
  mae = map_dbl(models, MAE,  cpts, data) 
  smae = map_dbl(models, spectral.error,  cpts, data) 
  result = tibble(
    method=c("Ratio", "Wang", "Galeano"),
    TDR=TDR, FDR=FDR, m_hat=m_hat, 
    MAE=mae, SMAE=smae, n=n, p=p, power=diff)
  #result = result %>% gather("metric", "value",-method,-n,-p,-power)
  #result = result %>% gather("metric", "value")
  return(result)
}


diff = 1 #this measures the strength of the jump
cp_num = 5 #the number of change-points in each cycle
num_repeats = 10 #how many periods we consider (note that the cyclic data is such that the last cycle is complete, not interrupted, can change it if necessary ) 
n = 100*cp_num #the length of each cycle
p = 3 #dimension of the time series (SPD dimension pxp)

#note that the total number of change-points is cp_num*num_repeats + num_repeats - 1
#because the end of the cycles are also change-points, n, 2n, 3n, ..., (num_repeats - 1)n

safe_run = safely(run.sim)
result=map(seq(1,1), ~safe_run(run+.x, n, p, cp_num, num_repeats, diff))
result

#m_hat is the absolute difference between the number of true and estimated changepoints
