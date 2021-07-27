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
#' @param age TRUE if covariate age should be included. simulated from N (40, 5)
#' @param sex TRUE if covariate sex should be included. simulated from Bernoulli (p), where p is a random number between 0 and 1.
#' @param n_homo number of homogeneous coefficients.
#' @param n_heter number of heterogeneous coefficients.
#' @param n_cont the number of continuous predictors. 0 if there will be no continuous predictors
#' @param n_SNP the number of SNP predictors. 0 if there will be no SNP predictors
#' @param rho a number between 0 and 1. controlling the degree of correlation between predictors
#' @param sigma a positive number. controlling the added noise to the simulated response variable
#' @param nlower the lower bound of the K sample sizes
#' @param nupper the upper bound of the K sample sizes
#' @param beta if the users wish to supply the coefficients on their own,
#' enter a coefficient matrix where the columns
#' containing the K coefficients of each predictor, for all genetic and non-genetic predictors.
#' @importFrom mvtnorm rmvnorm
#' @return the simulated data
#' \item{data}{the simulated data}
#' \item{beta}{the simulated beta matrix}
#' \item{homo_index}{the column numbers of homogeneous coefficients}
#' \item{heter_index}{the column numbers of heterogeneous coefficients}
#' @export
HomUHet.sim<-function(age=TRUE, sex=TRUE, n_homo, n_heter, K, n_cont, n_SNP,
                      rho=0.5, sigma=2,
                      nlower=50, nupper=300, beta=NULL) {




       if (n_SNP > 0) {

         allele_freq=runif(n_SNP, 0.05,0.5)
      }


      J=sum(c(isTRUE(age), isTRUE(sex),n_cont, n_SNP))
      t=0
      x=rep(0,J)
      y=0
      age=0
      sex=0
      study_label=0
      b=rep(0,40)
      n=0


      if(is.null(beta)){



      beta_gene=HomUHet.sim.beta(J=J-sum(c(isTRUE(age), isTRUE(sex))),K=K,n_homo=n_homo,
                            n_heter=n_heter, level=sample(c("l", "m", "h"),1))
      beta=beta_gene$beta
      homo_index=beta_gene$homo_index+sum(c(isTRUE(age), isTRUE(sex)))
      heter_index=beta_gene$heter_index+sum(c(isTRUE(age), isTRUE(sex)))

      ###### re organizing positions of homo heter when both cont and SNP preds exist

      if (n_cont > 0 & n_SNP > 0){
        beta_gene_re=matrix(0, ncol(beta_gene$beta), nrow=nrow(beta_gene$beta))

        homo_index_re=sample(seq(1:ncol(beta_gene_re)), n_homo)
        heter_index_re=sample(seq(1:ncol(beta_gene_re))[-homo_index_re], n_heter)
        beta_gene_re[,homo_index_re]=(beta_gene$beta)[,beta_gene$homo_index]
        beta_gene_re[,heter_index_re]=(beta_gene$beta)[,beta_gene$heter_index]
        beta=beta_gene_re
        homo_index=homo_index_re+sum(c(isTRUE(age), isTRUE(sex)))
        heter_index=heter_index_re+sum(c(isTRUE(age), isTRUE(sex)))
      }


      ####### binding age sex coefs

      beta=cbind(matrix(rep(rep(1,K), sum(c(isTRUE(age), isTRUE(sex)))), nrow=K),
                 beta)

      } else {



        check_func<-function(x){


          if(length(unique(x)) > 1){
            return(d="heter")
          } else {
            d=ifelse(unique(x)==0, "noise", "homo")

            return(d)
          }

        }

        check_index=apply(beta, 2, check_func)

        homo_index=which(check_index=="homo")
        heter_index=which(check_index=="heter")
      }

        while (t<K){

          n.temp=sample(seq(nlower,nupper),1)
          temp.x=matrix(0, ncol=J, nrow=n.temp)

          if (age==TRUE){
            temp_age=rnorm(n.temp, 40, 5)

          } else {
            temp_age=NULL
          }

          if (sex==TRUE) {
            temp_sex=rbinom(n.temp, 1, runif(1,0,1))

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
              divider=sqrt(allele_freq*(1-allele_freq))
              temp_snp_z1=rnorm(n.temp, 0, 1)
              temp_snp_z2=rnorm(n.temp, 0, 1)
              temp_snp_x1=ifelse(temp_snp_z1<=qnorm(allele_freq),1,0)
              temp_snp_x2=ifelse(temp_snp_z2<=qnorm(allele_freq),1,0)
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

              temp_snp_x1[,j]=ifelse(temp_snp_z1[,j]<=qnorm(allele_freq[j]),1,0)
              temp_snp_x2[,j]=ifelse(temp_snp_z2[,j]<=qnorm(allele_freq[j]),1,0)
              divider[j]=sqrt(allele_freq[j]*(1-allele_freq[j]))


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

      if (age==TRUE){
        age_name="age"
      } else {
        age_name=NULL
      }

      if (sex==TRUE){
        sex_name="sex"
      } else {
        sex_name=NULL
      }

      genetic_pred_num = J-sum(c(!is.null(age_name),!is.null(sex_name)))
      pred_names=c(age_name, sex_name,paste(rep("Pred_V",genetic_pred_num), seq(1:genetic_pred_num),sep="" ))

      names(data)[-(1:2)]=pred_names

       #####
      homo_index=homo_index+2
      heter_index=heter_index+2

      list(data=data, beta=beta, homo_index=homo_index,
           heter_index=heter_index)


    }




#' simulates homogeneous and heterogeneous coefficients of predictors
#'
#' this function outputs matrix of coefficients of predictors
#'
#' @param J The total number of predictors including the predictors with homogeneous and heterogeneous effects and the predictors without effects
#' @param K The number of studies
#' @param n_homo the number of homogeneous coefficients
#' @param n_heter the number of heterogeneous coefficients
#' @param level the level of heterogeneity in the heterogeneous coefficients,
#' where "l" stands for low, "m" stands for medium and "h" stands for high.
#'
#' @return the simulated coefficient matrix and miscellenance information about it
#'  \item{beta}{the K x J coefficient matrix}
#'  \item{J}{the number of predictors including both predictors which have effects and which do not}
#'  \item{K}{the number of studies}
#'  \item{homo_index}{a vector containing the column numbers of homogeneous coefficients in the coefficient matrix}
#'  \item{heter_index}{a vector containing the column numbers of homogeneous coefficients in the coefficient matrix}
#' @export
HomUHet.sim.beta<-function(J, K, n_homo, n_heter, level=c("l","m","h")){

  # divide homo:

  n_homo_1=floor(n_homo/2)

  n_homo_2=n_homo-floor(n_homo/2)


  # divide hetero cluster:
  cluster_1=floor(n_heter/3)
  cluster_2=floor(n_heter/3)
  cluster_3=n_heter-2*(floor(n_heter/3))

  # heter, divide K:

  # cluster 1

  n_heter_c1_1=floor(K/2)
  n_heter_c1_2=K-floor(K/2)

  # cluster 2 and 3 :
  n_heter_c2_1=floor(K/3)
  n_heter_c2_2=floor(K/3)
  n_heter_c2_3=K-2*floor(K/3)



  l_param=c(1,3,-3,-1,
            1,1.5,-1.5,-1,
            1,1.5,2,2.5,3,3.5,
            -1.5,-1,-2.5,-2,-3.5,-3)
  m_param=c(1,3,-3,-1,
            2,2.5,-2.5,-2,
            1,1.5,3,3.5,5,5.5,
            -1.5,-1,-3.5,-3,-5.5,-5)

  h_param=c(3,6,-6,-3,
            5,6.5,-6.5,-5,
            1,2,5,6,9,10,
            -2,-1,-6,-5,-10,-9)

  if (level=="l"){
    param=l_param
  } else if (level=="m") {
    param=m_param
  } else if (level=="h") {
    param=h_param
  }


  homo=c(runif(n_homo_1,param[1],param[2]),runif(n_homo_2,param[3],param[4]))

  homo=t(matrix(rep(homo,K), ncol=K))

  heter=matrix(0, ncol=n_heter, nrow=K)


  for (j in 1:cluster_1){
    heter[,j]=sample(c(runif(n_heter_c1_1, param[5],param[6]), runif(n_heter_c1_2, param[7],param[8])),K)
  }

  for (j in (cluster_1+1):(2*cluster_1)){
    heter[,j]=sample(c(runif(n_heter_c2_1, param[9],param[10]), runif(n_heter_c2_2, param[11],param[12]), runif(n_heter_c2_3, param[13],param[14])),K)
  }

  for (j in (2*cluster_1+1):(2*cluster_1+cluster_3)){
    heter[,j]=sample(c(runif(n_heter_c2_1, param[15],param[16]), runif(n_heter_c2_2, param[17],param[18]), runif(n_heter_c2_3, param[19],param[20])),K)
  }

  number_col=n_homo_1+floor((1/14)*J)+n_homo_2+floor((1/14)*J)+2*n_heter

  beta=matrix(0,ncol=J, nrow=K)

  if (n_homo > 0){
    homo_index=c(1:n_homo_1,c((n_homo_1+floor((1/14)*J)+1):(n_homo_1+floor((1/14)*J)+n_homo_2)))
  } else {
    homo_index=NULL
  }

  if (n_heter > 0){
    heter_index=seq(n_homo+2*floor((1/14)*J)+1, n_homo+2*floor((1/14)*J)+n_heter*2, by=2)
  } else {
    heter_index=NULL
  }


  beta[,homo_index]=homo

  beta[,heter_index]=heter



  list(J=J, K=K, beta=beta,
       homo_index=homo_index, heter_index=heter_index)


}
