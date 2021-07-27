#'
#' fit a two-step penalized regression model
#'
#' This function outputs the names of predictors with homogeneous or heterogeneous predictors across multiple data sets, the estimates of predictors, and solution plots
#'
#' @param data The input dataframe containing observations from all studies where the first column is the study label,
#' the second column is the response variable and the following columns are the predictors.
#' @param solution_path_plot TRUE if outputting solution path plots is desired
#' @return the names of identified predictors and their estimated effects
#'  \item{Homo}{a character string of names of homogeneous predictors}
#'  \item{Heter}{a character string of N=names of heterogeneous predictors}
#'  \item{coefficients}{a data frame containing estimated coefficients of the homogeneous and heterogeneous predictors in K studies}
#'
#' @importFrom dplyr arrange group_by summarise_all
#' @importFrom graphics abline plot
#' @importFrom stats coef lm logLik var
#' @importFrom mvtnorm rmvnorm
#' @import glmnet
#' @import gglasso
#'
#'
#' @export
HomUHet<-function(data,solution_path_plot=FALSE){
  ##### sorting data by study

  if (is.null(names(data))){
    pred_names=past("V",c(3:ncol(data)), sep="")
  }
  pred_names=names(data)[-c(1:2)]
  x=as.matrix(data[,-c(1:2)])
  y=data[,2]
  study_label=data[,1]
  n=as.data.frame(table(study_label))[,2]
  data=cbind(y,x)
  data=as.data.frame(cbind(study_label,data))
  data=dplyr::arrange(data,study_label)

  J=ncol(x)
  K=length(unique(study_label))
  lambda=c(0,0)
  lambda[1]=ifelse(sum(n)>J,0.0001,0.01)
  lambda[2]=ifelse(sum(n)>(J*K),0.001,0.05)






  #### fitting adaptive LASSO with BIC

  ini.beta<-glmnet::cv.glmnet(as.matrix(x),y,alpha=0,standardize=FALSE)

  beta.bic.fit<-glmnet::glmnet(as.matrix(x),y,alpha=1,penalty.factor=1/abs(coef(ini.beta,s="lambda.min")[-1]),standardize=FALSE,lambda.min.ratio = lambda[1])
  Bic<-function(y,x,beta,df){
    bic=length(y)*log(sum((y-x%*%beta)^2)/length(y))+log(length(y))*df
  }
  bic_score=rep(0,length(beta.bic.fit$df))
  for (q in 1:length(beta.bic.fit$df)){
    bic_score[q]=Bic(y,x,beta.bic.fit$beta[,q],beta.bic.fit$df[q])
  }
  Bic_optimal=which.min(bic_score)
  non.zero.beta=which(beta.bic.fit$beta[,Bic_optimal]!=0)
  zero.beta=which(beta.bic.fit$beta[,Bic_optimal]==0)


  #### finding OLS estimates for beta j
  data=as.data.frame(cbind(y,x[,non.zero.beta]))
  beta.ls.fit<-lm(y~.,data=data)

  y.res=y-cbind(rep(1,length(y)),x[,non.zero.beta])%*%beta.ls.fit$coefficients

  ### step 1 estimates
  step1_estimates=beta.bic.fit$beta[,Bic_optimal]

  step1_estimates[non.zero.beta]=beta.ls.fit$coefficients[-1]

  ##### creating joint data set


  x.by.var<-matrix(0,nrow=sum(n),ncol=J*K)

  for (j in 1:J){

    for (i in 1:K){
      temp=x.by.var[,((j-1)*K+1):(j*K)]

      if(i<2){
        temp[1:n[i],i]=x[1:n[i],j]
      }

      else
        temp[(sum(n[seq(1:(i-1))])+1):(sum(n[seq(1:(i))])),i]=x[(sum(n[seq(1:(i-1))])+1):(sum(n[seq(1:(i))])),j]

      x.by.var[,((j-1)*K+1):(j*K)]=temp
    }
  }

  ##### fitting group lasso



  group=rep(0,J*K)
  for (j in 1:J){

    group[((j-1)*K+1):(j*K)]=rep(j,K)

  }



  #### obtaining preliminary fit for delta j

  ini.delta.cv.fit<-glmnet::cv.glmnet(as.matrix(x.by.var),y.res,alpha=0,standardize=FALSE,intercept=TRUE)
  ini.delta.fit<-glmnet::glmnet(as.matrix(x.by.var),y.res,alpha=0,standardize=FALSE,intercept=TRUE,
                        lambda=ini.delta.cv.fit$lambda.min)

  ini.delta.norm=rep(0,J)
  ini.delta=ini.delta.fit$beta
  for (j in 1:J){
    ini.delta.norm[j]=norm(ini.delta[((j-1)*K+1):(j*K)],type="2")

  }




  #### fitting with all delta-j






  ebic.grp.fit<-gglasso::gglasso (as.matrix(x.by.var), y.res, group = group,
                         loss = "ls",
                         nlambda = 100,
                         lambda.factor = lambda[2],
                         pf = rep(sqrt(K),J)/(ini.delta.norm+0.000000000001),
                         dfmax = 1400,
                         pmax = 1400,
                         eps = 1e-08, maxit = 3e+08, intercept=TRUE)



  ####creating separate dataset

  data.by.study<-list()

  data.by.study[[1]]=cbind(y.res[1:(n[1])],x[1:(n[1]),])


  for (k in 2:K){

    data.by.study[[k]]=cbind(y.res[(sum(n[1:(k-1)])+1):sum(n[1:k])],x[(sum(n[1:(k-1)])+1):sum(n[1:k]),])
  }




  third.term.grp<-function(x,J,K,gamma){

    num.p=length(x)/K

    if (num.p==0){
      factorial.part=0
    }
    else
      factorial.part=sum(log(seq(1:J)))-(sum(log(seq(1:num.p)))+sum(log(seq(1:(J-num.p)))))
    third.term=2*gamma*factorial.part
    return(third.term)
  }


  index=abs(ebic.grp.fit$beta)>0
  L=which.max(which(ebic.grp.fit$df<((n[which.min(n)]*2/3)*K)))
  ebic=rep(0,L)

  for (l in 1:length(ebic)){

    temp.select_index=unique(group[index[,l]])

    temp.fit=as.data.frame(cbind(as.factor(group),(ebic.grp.fit$beta)[,l]))

    temp.selection=temp.fit[temp.fit$V1%in%temp.select_index,2]

    temp.third=third.term.grp(temp.selection,J,K,1)

    temp.ls.fit.lik=rep(0,K)

    if (l<2){

      for (k in 1:K){
        temp.data=as.data.frame((data.by.study[[k]])[,1])
        names(temp.data)=c("V1")
        temp.ls.fit=lm(V1~.,data=temp.data)
        temp.ls.fit.lik[k]=logLik(temp.ls.fit)
      }


    }else
    {
      temp.ls.coef=matrix(0,nrow=K,ncol=length(temp.select_index))
      for (k in 1:K){
        temp.data=as.data.frame(cbind((data.by.study[[k]])[,1],((data.by.study[[k]])[,-1])[,temp.select_index]))
        temp.ls.fit=lm(V1~.,data=temp.data)
        temp.ls.fit.lik[k]=logLik(temp.ls.fit)
        temp.ls.coef[k,]=temp.ls.fit$coefficients[-1]
      }

    }


    temp.dev=-2*(sum(temp.ls.fit.lik))
    norm.lasso=rep(0,length(temp.selection)/K)
    norm.ls=rep(0,length(temp.selection)/K)
    if (length(temp.selection)/K>0){
      for (i in 1:(length(temp.selection)/K)){

        norm.lasso[i]=norm(temp.selection[((i-1)*K+1):(i*K)],type="2")

        temp.ls.coef[is.na(temp.ls.coef[,i]),i]=0### replace NA by 0, ls fit model saturation
        norm.ls[i]=norm(temp.ls.coef[,i],type="2")
      }
      temp.df=length(temp.selection)/K+(K-1)*sum(norm.lasso/norm.ls)
    }else

    {temp.df=0}

    ebic[l]=temp.dev+temp.df*log(length(y.res))+temp.third

  }






  ebic.lambda=ebic.grp.fit$lambda[which.min(ebic)]

  ebic.coef=ebic.grp.fit$beta[,which.min(ebic)]

  ebic.fitted.delta=as.data.frame(cbind(as.factor(group),ebic.coef))


  ebic.selection=ebic.fitted.delta%>%dplyr::group_by(V1)%>%dplyr::summarise_all(var)### sd not working




  #### assessing performance

  Homo=intersect(non.zero.beta,which(ebic.selection[,2]==0))
  Homo_names=pred_names[Homo]
  Heter=which(ebic.selection[,2]!=0)
  Heter_names=pred_names[Heter]
  noise=intersect(zero.beta,which(ebic.selection[,2]==0))

  #### estimates
  if (length(Homo)>0){
    Homo_estimates=matrix(0,ncol=(K+1),nrow=length(Homo))
    for (r in 1:length(Homo)){
      Homo_estimates[r,]=c(Homo_names[r],rep(step1_estimates[Homo[r]],K))
    }
    Homo_estimates=as.matrix(Homo_estimates)
    Homo_estimates=cbind(Homo_estimates,rep("Homogeneous", nrow(Homo_estimates)))
  } else {
    Homo_estimates=NULL
  }

  if (length(Heter)>0){
    Heter_estimates=matrix(0,ncol=(K+1),nrow=length(Heter))

    for (p in 1:length(Heter)){
      Heter_estimates[p,]=c(Heter_names[p],ebic.fitted.delta[ebic.fitted.delta$V1==Heter[p],2]+step1_estimates[Heter[p]])
    }

    Heter_estimates=as.matrix(Heter_estimates)
    Heter_estimates=cbind(Heter_estimates,rep("Heterogeneous", nrow(Heter_estimates)))
  } else {
    Heter_estimates=NULL
  }


  coefficients=rbind(Homo_estimates,
                       Heter_estimates)
  coefficients=as.data.frame(coefficients)
  names(coefficients)=c("predictor",paste("study",seq(1:K),sep=""),"type")

if (solution_path_plot==TRUE){

  y_name=names(y)
  plot(beta.bic.fit,xvar = "lambda",main=paste(y_name,"step1",sep = " "))
  abline(v=log(beta.bic.fit$lambda[Bic_optimal]))
  plot(ebic.grp.fit,main=paste(y_name,"step2",sep = " "))
  abline(v =log(ebic.lambda))
}



  list(Homo=Homo_names,
       Heter=Heter_names,
       coefficients=coefficients)
}

#' simulate multiple data sets with both homogeneous and heterogeneous effects from the predictors
#'
#' this function simulate data
#'
#' @param J the total number of predictors.
#' @param K the number of studies.
#' @param param_age should be a vector containing the mean of age, followed by the standard deviation.
#' @param param_sex should be a number between 0 and 1 indicating the proportion of one sex
#' @param n_cont the number of continuous predictors. 0 if there will be no continuous predictors
#' @param n_SNP the number of SNP predictors. 0 if there will be no SNP predictors
#' @param beta the K x J coefficient matrix
#' @param rho a number between 0 and 1. controlling the degree of correlation between predictors
#' @param sigma a positive number. controlling the added noise to the simulated response variable
#' @param allele_freq a vector of length n_SNP containing the allele frequencies for the SNPs.
#' @param nlower the lower bound of the K sample sizes
#' @param nupper the upper bound of the K sample sizes
#' @importFrom mvtnorm rmvnorm
#' @return the simulated data
#' @export
HomUHet.sim<-function(J, K, param_age=NULL, param_sex=NULL, n_cont, n_SNP, allele_freq, beta,
                      rho=0.5,sigma=2,
                      nlower=50,nupper=300 ) {


   if (is.null(beta)){
      stop("beta is missing")
    } else {

      if(K != nrow(beta) ){
        stop ("K does not equal to the number of rows of beta")
      } else if (J != ncol(beta)) {
        stop ("J does not equal to the number of columns of beta")
      } else if (n_SNP != length(allele_freq)) {
        stop ("n_SNP needs to match the length of allele_freq")
      }

      K=nrow(beta)
      J=ncol(beta)
      t=0
      x=rep(0,J)
      y=0
      age=0
      sex=0
      study_label=0
      b=rep(0,40)
      n=0

      #J, K, param_age=c(40, 5), param_sex=0.5, n_cont, n_SNP


        while (t<K){

          n.temp=sample(seq(nlower,nupper),1)
          temp.x=matrix(0, ncol=J, nrow=n.temp)

          if (!is.null(param_age)){
            temp_age=rnorm(n.temp, param_age[1], param_age[2])

          } else {
            temp_age=NULL
          }

          if (!is.null(param_sex)) {
            temp_sex=rbinom(temp.n, 1, param_sex)

          } else {
            temp_sex= NULL
          }


          # Cont

          if (n_cont > 0){

            if (n_cont==1){
              temp_cont_x=rnorm(n.temp, 0, 1)
            }

            v=c(1,rep(NA,((n_cont)-1)))

            for (i in 2:(n_cont)){
              v[i]=rho^(i-1)
            }

            Sigma=toeplitz(v)

            temp_cont_x=mvtnorm::rmvnorm(n.temp, mean = rep(0, n_cont), sigma = Sigma,
                                    method="svd", pre0.9_9994 = FALSE, checkSymmetry = TRUE)

          } else {
            temp_cont_x=NULL
          }

          # SNP

          if (n_SNP > 0) {

            if (n_SNP==1){
              divider=sqrt(allele.freq*(1-allele.freq))
              temp_snp_z1=rnorm(n.temp, 0, 1)
              temp_snp_z2=rnorm(n.temp, 0, 1)
              temp_snp_x1=ifelse(temp_snp_z1<=qnorm(allele.freq),1,0)
              temp_snp_x2=ifelse(temp_snp_z2<=qnorm(allele.freq),1,0)
              temp_snp_x=(temp_snp_x1+temp_snp_x2)*divider
            }

            v=c(1,rep(NA,((n_SNP)-1)))

            for (i in 2:(n_SNP)){
              v[i]=rho^(i-1)
            }

            Sigma=toeplitz(v)

            divider=rep(0,n_SNP)
            temp_snp_x1=matrix(0,ncol=n_SNP,nrow=n.temp)
            temp_snp_x2=matrix(0,ncol=n_SNP,nrow=n.temp)


            temp_snp_z1=mvtnorm::rmvnorm(n.temp, mean = rep(0, n_SNP), sigma = Sigma,
                                     method="svd", pre0.9_9994 = FALSE, checkSymmetry = TRUE)
            temp_snp_z2=mvtnorm::rmvnorm(n.temp, mean = rep(0, n_SNP), sigma = Sigma,
                                     method="svd", pre0.9_9994 = FALSE, checkSymmetry = TRUE)

            for (j in 1:n_SNP){

              temp_snp_x1[,j]=ifelse(temp_snp_z1[,j]<=qnorm(allele.freq[j]),1,0)
              temp_snp_x2[,j]=ifelse(temp_snp_z2[,j]<=qnorm(allele.freq[j]),1,0)
              divider[j]=sqrt(allele.freq[j]*(1-allele.freq[j]))


            }

            temp_snp_x=(temp_snp_x1+temp_snp_x2)*divider

          } else {

            temp_snp_x=NULL
          }


          temp.x=cbind(temp_age, temp_sex, temp_cont_x, temp_snp_x)
          temp.x=cbind(temp.x, matrix(0, ncol=J-ncol(temp.x), nrow=n.temp))

          study_label.temp=rep(t+1,n.temp)


          temp.rb=beta[(t+1),]
          temp.y=rnorm(1,0,5)+temp.x%*%temp.rb+rnorm(n.temp,0,sigma)

          x=rbind(x,scale(temp.x))
          y=c(y,temp.y)
          n=c(n,n.temp)
          study_label=c(study_label,study_label.temp)
          t=t+1
        }

      x=x[-1,]
      y=y[-1]
      n=n[-1]
      study_label=study_label[-1]

      data=as.data.frame(cbind(study_label, y, x))

      if (! is.null(param_age)){
        age_name="age"
      } else {
        age_name=NULL
      }

      if (! is.null(param_sex)) {
        sex_name="sex"
      } else {
        sex_name=NULL
      }

      pred_names=c(age_name, sex_name,paste(rep("Pred_V",ncol(x)), seq(1:ncol(x)),sep="" ))

      names(data)[-(1:2)]=pred_names

      list(data)

        }



    }




#' simulates homogeneous and heterogeneous coefficients of predictors
#'
#' this function outputs matrix of coefficients of predictors
#'
#' @param J The total number of predictors including the predictors with homogeneous and heterogeneous effects and the predictors without effects
#' @param K The number of studies
#' @param n_homo the number of homogeneous coefficients
#' @param n_heter the number of heterogeneous coefficients
#' @param param_homo the two parameters of the uniform distribution used to simulate the homogeneous coefficients
#' @param param_heter the two parameters of the uniform distribution used to simulate the heterogeneous coefficients
#' @param example provides some examples of coefficient matrix for various number of studies (4 or 10) and levels of heterogeneity.
#' enter the level of heterogeneity in the first entry. "l", "m", and "h" represent the level of heterogeneity in coefficients
#' enter the number of studies in the second entry. 4 or 10.
#' where "l" stands for low, "m" stands for medium and "h" stands for high.
#'
#' @return the simulated coefficient matrix and miscellenance information about it
#'  \item{beta}{the K x J coefficient matrix}
#'  \item{J}{the number of predictors including both predictors which have effects and which do not}
#'  \item{K}{the number of studies}
#'  \item{homo_index}{a vector containing the column numbers of homogeneous coefficients in the coefficient matrix}
#'  \item{heter_index}{a vector containing the column numbers of homogeneous coefficients in the coefficient matrix}
#' @export
HomUHet.sim.beta<-function(J, K, n_homo, n_heter, param_homo=c(-10,10),
                           param_heter=c(-10,10), example=c(c("l","m","h"),c(4,10))){

  #if (J==1400 & K %in% c(4,10) & n_homo=10 &n_heter=30 &is.null(param_homo)
  #    & is.null(param_heter) & !is.null(param_default) ){

  if (!is.null(example)){
  J=1400
  K=example[2]
  level=example[1]

  if (K==4){
    if (level=="h"){
      fixed.b=c(runif(5,1,3),runif(5,-3,-1))
      b.matrix=matrix(0,ncol=30,nrow=K)

      for (j in 1:30){
        if (j<11){
          b.matrix[,j]=sample(c(runif(2,1,1.5),runif(2,-1.5,-1)),K)
        }

        else {if (j<21){
          b.matrix[,j]=sample(c(runif(1,1,1.5),runif(1,2,2.5),runif(2,3,3.5)),K)
        }
          else {
            b.matrix[,j]=sample(c(runif(1,-1.5,-1),runif(1,-2.5,-2),runif(2,-3.5,-3)),K)
          }
        }
      }
      modified.b.matrix<-matrix(0,ncol=60,nrow=K)

      for (j in 1:30){

        modified.b.matrix[,(j*2-1)]=b.matrix[,j]
      }

    } else if (level=="m"){
      fixed.b=c(runif(5,1,3),runif(5,-3,-1))
      b.matrix=matrix(0,ncol=30,nrow=K)

      for (j in 1:30){
        if (j<11){
          b.matrix[,j]=sample(c(runif(2,2,2.5),runif(2,-2.5,-2)),K)
        }

        else {if (j<21){
          b.matrix[,j]=sample(c(runif(1,1,1.5),runif(1,3,3.5),runif(2,5,5.5)),K)
        }
          else {
            b.matrix[,j]=sample(c(runif(1,-1.5,-1),runif(1,-3.5,-3),runif(2,-5.5,-5)),K)
          }
        }
      }

      modified.b.matrix<-matrix(0,ncol=60,nrow=K)

      for (j in 1:30){

        modified.b.matrix[,(j*2-1)]=b.matrix[,j]
      }

    } else {
      fixed.b=c(runif(5,3,6),runif(5,-6,-3))
      b.matrix=matrix(0,ncol=30,nrow=K)

      for (j in 1:30){
        if (j<11){
          b.matrix[,j]=sample(c(runif(2,5,6.5),runif(2,-6.5,-5)),K)
        }

        else {if (j<21){
          b.matrix[,j]=sample(c(runif(1,1,2),runif(1,5,6),runif(2,9,10)),K)
        }
          else {
            b.matrix[,j]=sample(c(runif(1,-2,-1),runif(1,-6,-5),runif(2,-10,-9)),K)
          }
        }
      }
      modified.b.matrix<-matrix(0,ncol=60,nrow=K)

      for (j in 1:30){

        modified.b.matrix[,(j*2-1)]=b.matrix[,j]
      }

    }} else {

      if (level=="h"){
        fixed.b=c(runif(5,1,3),runif(5,-3,-1))
        b.matrix=matrix(0,ncol=30,nrow=K)

        for (j in 1:30){
          if (j<11){
            b.matrix[,j]=sample(c(runif(5,1,1.5),runif(5,-1.5,-1)),K)
          }

          else {if (j<21){
            b.matrix[,j]=sample(c(runif(3,1,1.5),runif(3,2,2.5),runif(4,3,3.5)),K)
          }
            else {
              b.matrix[,j]=sample(c(runif(3,-1.5,-1),runif(3,-2.5,-2),runif(4,-3.5,-3)),K)
            }
          }
        }
        modified.b.matrix<-matrix(0,ncol=60,nrow=K)

        for (j in 1:30){

          modified.b.matrix[,(j*2-1)]=b.matrix[,j]
        }
      } else if (level=="m"){
        fixed.b=c(runif(5,1,3),runif(5,-3,-1))
        b.matrix=matrix(0,ncol=30,nrow=K)

        for (j in 1:30){
          if (j<11){
            b.matrix[,j]=sample(c(runif(5,2,2.5),runif(5,-2.5,-2)),K)
          }

          else {if (j<21){
            b.matrix[,j]=sample(c(runif(3,1,1.5),runif(3,3,3.5),runif(4,5,5.5)),K)
          }
            else {
              b.matrix[,j]=sample(c(runif(3,-1.5,-1),runif(3,-3.5,-3),runif(4,-5.5,-5)),K)
            }
          }
        }

        modified.b.matrix<-matrix(0,ncol=60,nrow=K)

        for (j in 1:30){

          modified.b.matrix[,(j*2-1)]=b.matrix[,j]
        }
      } else {
        fixed.b=c(runif(5,1,6),runif(5,-6,-1))
        b.matrix=matrix(0,ncol=30,nrow=K)

        for (j in 1:30){
          if (j<11){
            b.matrix[,j]=sample(c(runif(5,5,6.5),runif(5,-6.5,-5)),K)
          }

          else {if (j<21){
            b.matrix[,j]=sample(c(runif(3,1,2),runif(3,5,6),runif(4,9,10)),K)
          }
            else {
              b.matrix[,j]=sample(c(runif(3,-2,-1),runif(3,-6,-5),runif(4,-10,-9)),K)
            }
          }
        }
        modified.b.matrix<-matrix(0,ncol=60,nrow=K)

        for (j in 1:30){

          modified.b.matrix[,(j*2-1)]=b.matrix[,j]
        }
      }

    }
    #c(fixed.b[1:5],rep(0,100),fixed.b[6:10],rep(0,100),modified.b.matrix[(t+1),])
    beta=cbind(t(matrix(rep(fixed.b[1:5],K), ncol=K)), matrix(0,nrow=K, ncol=100),
          t(matrix(rep(fixed.b[6:10],K), ncol=K)), matrix(0,nrow=K, ncol=100),
          modified.b.matrix, matrix(0,nrow=K, ncol=J-270))

    homo_index=c(c(1:5), c(106:110))
    heter_index=210+seq(1,60,by=2)

  } else {



    if (J < 1) {
    stop ("number of predictors is less than 1")

  }

  if (K <= 1) {
    stop ("number of studies is less than or equal to 1")
  }

  if (sum(c(n_homo,n_heter)) > J) {
    stop ("number of predictors is less than the total number of homogeneous and
          heterogeneous predictors")
  }


  if (param_homo[1] > param_homo[2]){

    stop ("incorrect lower and upper bound for simulating homogeneous effects")
  }

  if (param_heter[1] >= param_heter[2]){

    stop (
      "incorrect lower and upper bound for simulating heterogeneous
     effects. upper bound should be more than lower bound"
    )
  }

  beta=matrix(0, ncol=J, nrow=K)

  if (n_homo > 0 & n_heter > 0) {
    homo_values = runif(n_homo, param_homo[1], param_homo[2])
    beta_homo = t(matrix(rep(homo_values, K), ncol=K))
    homo_index = sample(seq(1:J), n_homo)
    beta[,homo_index]=beta_homo

      heter_gen<-function(x){

        return(runif(x[1], min=x[2], max=x[3]))
      }

      heter_coef_param=rbind(rep(K, n_heter), matrix(rep(param_heter, n_heter), ncol= n_heter))
      beta_heter=apply(heter_coef_param, 2, heter_gen)

      heter_index= sample ((seq(1:J))[-c(homo_index)], n_heter)
      beta[,heter_index]=beta_heter


    } else if (n_homo <= 0 & n_heter > 0) {
      heter_gen<-function(x){

        return(runif(x[1], min=x[2], max=x[3]))
      }

      heter_coef_param=rbind(rep(K, n_heter), matrix(rep(param_heter, n_heter), ncol= n_heter))
      beta_heter=apply(heter_coef_param, 2, heter_gen)
      heter_index= sample ((seq(1:J)), n_heter)
      beta[,heter_index]=beta_heter
    } else if (n_homo > 0 & n_heter <= 0){
      homo_values = runif(n_homo, param_homo[1], param_homo[2])
      beta_homo = t(matrix(rep(homo_values, K), ncol=K))
      homo_index = sample(seq(1:J), n_homo)
      beta[,homo_index]=beta_homo
    } else {

      stop ("no information about homogenous or heterogenous effects are supplied")
    }
  }



  #beta=cbind(beta_homo, beta_heter, matrix(0, ncol=(J-sum(c(n_homo,n_heter))), nrow = K)


  list(J=J, K=K, beta=beta,
       homo_index=homo_index, heter_index=heter_index)


}
