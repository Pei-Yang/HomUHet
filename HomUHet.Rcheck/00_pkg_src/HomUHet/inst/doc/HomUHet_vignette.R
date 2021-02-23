## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include=FALSE-----------------------------------------------------------
options(tinytex.verbose = TRUE)

## -----------------------------------------------------------------------------
library(HomUHet)
data(HomUHet_data)
HomUHet_data[1:5,1:4]
# input data of 4 studies, same set of 500 predictors for each study. 
#The first column is the study label, the second column is the response variable y, 
#and Pred_V1 and Pred_V2 are the first 2 predictors.

unique(HomUHet_data$study_label) 
# study_label contains 4 different integers labeling the studies. 

x1=HomUHet_data[HomUHet_data$study_label==1,] 
# extracting observations from study 1
(colMeans(x1[,-(1:2)]))[1:5] # checking if the means are 0. showing the first 5 predictors here.
(apply(x1[,-(1:2)],2,sd))[1:5] # checking if the variances are 1. showing first 5 here.


## -----------------------------------------------------------------------------
fit=HomUHet(data=HomUHet_data, solution_path_plot = TRUE)
fit$Homo # the identified homogeneous predictors (names are not provided so these are column names)
fit$Heter # the identified heterogeneous predictors
(fit$coefficients[(fit$coefficients)$type=="Homogeneous",])[1:5,] 
# some of the estimated homogeneous coefficients 
(fit$coefficients[(fit$coefficients)$type=="Heterogeneous",])[1:5,] 
# some of the estimated coefficients


## -----------------------------------------------------------------------------
beta=cbind(c(2,2,2), matrix(0,ncol=10,nrow=3),c(2,4,6))
# constructing a beta matrix
K=nrow(beta) # the number of studies
J=ncol(beta) # the number of predictors
allele_freq=runif(J,0.05,0.5) # allele frequency 
mydata=HomUHet.sim(Pred_type="SNP", J=J, K=K, beta=beta,
                      rho=0.5,sigma=2,
                      nlower=50,nupper=200, allele_freq=allele_freq)

mydata[1:2,1:5] # showing the first few lines and columns of the output data.


