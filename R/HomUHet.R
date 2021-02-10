
#' HomUHet
#'
#' This function outputs the names of predictors with homogeneous or heterogeneous predictors across multiple data sets, the estimates of predictors, and solution plots
#'
#' @param data the data containing observations from all studies, with the response variable in the first column, followed by the predictor matrix.
#' @param sid the study id for each observation in data
#' @param solution_path TRUE if the user wishes to output the solution path plots. Default is FALSE
#' @param y_name if needed, a response variable name for the solution path plots. Default is NULL.
#' @return names of Homogeneous and Heterogeneous predictors, estimates of predictors, solution path plots
#' 
#' @importFrom dplyr arrange group_by summarise_all
#' @importFrom HDeconometrics ic.glmnet
#' @importFrom graphics abline plot
#' @importFrom stats coef lm logLik var
#' @import glmnet
#' @import gglasso
#'
#' 
#' @export
HomUHet<-function(data,sid,solution_path=FALSE,y_name=NULL){
  ##### sorting data by study
  
  data=as.data.frame(cbind(sid,data))
  data=dplyr::arrange(data,sid)
  n=as.data.frame(table(sid))[,2]
  x=as.matrix(data[,-c(1:2)])
  y=data[,2]
  sid=data[,1]
  J=ncol(x)
  K=length(unique(sid))
  lambda=c(0,0)
  lambda[1]=ifelse(sum(n)>J,0.0001,0.01)
  lambda[2]=ifelse(sum(n)>(J*K),0.001,0.05)
  
  
  

  
 
  #### fitting adaptive LASSO with BIC

  ini.beta<-glmnet::cv.glmnet(as.matrix(x),y,alpha=0,standardize=FALSE)
  beta.bic.fit<-HDeconometrics::ic.glmnet(as.matrix(x),y,alpha=1,crit="bic",penalty.factor=1/abs(coef(ini.beta,s="lambda.min")[-1]),standardize=FALSE,lambda.min.ratio = lambda[1])
  non.zero.beta=which(beta.bic.fit$coefficients[-1]!=0)
  zero.beta=which(beta.bic.fit$coefficients[-1]==0)
  
  
  #### finding OLS estimates for beta j
  data=as.data.frame(cbind(y,x[,non.zero.beta]))
  beta.ls.fit<-lm(y~.,data=data)
  
  y.res=y-cbind(rep(1,length(y)),x[,non.zero.beta])%*%beta.ls.fit$coefficients
  
  ### step 1 estimates
  step1_estimates=beta.bic.fit$coefficients[-1]
  
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
  Heter=which(ebic.selection[,2]!=0)
  noise=intersect(zero.beta,which(ebic.selection[,2]==0))
  
  #### estimates
  if (length(Homo)>0){
    Homo_estimates=matrix(0,ncol=(K+1),nrow=length(Homo))
    for (r in 1:length(Homo)){
      Homo_estimates[r,]=c(Homo[r],rep(step1_estimates[Homo[r]],K))
    }
    Homo_estimates=as.matrix(Homo_estimates)
  } else {
    Homo_estimates=NULL
  }
  
  if (length(Heter)>0){
    Heter_estimates=matrix(0,ncol=(K+1),nrow=length(Heter))
    
    for (p in 1:length(Heter)){
      Heter_estimates[p,]=c(Heter[p],ebic.fitted.delta[ebic.fitted.delta$V1==Heter[p],2]+step1_estimates[Heter[p]])
    } 
    
    Heter_estimates=as.matrix(Heter_estimates)
  } else {
    Heter_estimates=NULL
  }
  
  
  #### outputting solution path
  
  if (solution_path==TRUE){
    plot(beta.bic.fit$glmnet,xvar = "lambda",main=paste(y_name,"step1",sep = " "))
    abline(v=log(beta.bic.fit$lambda))
    plot(ebic.grp.fit,main=paste(y_name,"step2",sep = " "))
    abline(v =log(ebic.lambda))
  }
  list(Homo=Homo,
       Heter=Heter,
       Homo_estimates=Homo_estimates,
       Heter_estimates=Heter_estimates,
       Solution_Path=Solution_Path)
}