## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------


## -----------------------------------------------------------------------------
library(HomUHet)
data(HomUHet_data)
HomUHet_data[1:5,1:4]# input data of 4 studies, same set of 500 predictors for each study.

dim(HomUHet_data)

unique(HomUHet_data$study_label) # study_label contains 4 different integers labeling the studies. 

x1=HomUHet_data[HomUHet_data$study_label==1,] # extracting observations from study 1
(colMeans(x1[,-(1:2)]))[1:5] # checking if the means are 0. showing the first 5 predictors here.
(apply(x1[,-(1:2)],2,sd))[1:5] # checking if the variances are 1. showing first 5 here.


## -----------------------------------------------------------------------------
 HomUHet(data=HomUHet_data,solution_path_plot=TRUE)

## -----------------------------------------------------------------------------
fit=HomUHet(data=HomUHet_data)
fit$Homo # the identified homogeneous predictors
fit$Heter # the identified heterogeneous predictors
fit$coefficients[1:2,] # the estimated coefficients
fit$coefficients[15:17,] # the estimated coefficients


