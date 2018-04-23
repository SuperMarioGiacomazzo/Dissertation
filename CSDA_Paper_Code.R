#########################################################
#Paper:"REGULARIZATION METHODS FOR SUBSET ARMA SELECTION"
#Authors:Mario Giacomazzo (Arizona State University) 
#        Yiannis Kamarianakis (Arizona State University)
#Year:2018
#########################################################


###############################################################################
#Required R Packages
###############################################################################
library(doParallel) #Needed for Potential Parallel Processing
library(foreach) #Needed for Potential Parallel Processing
library(glmnet)  #Needed for Adaptive Lasso/Adaptive Elastic Net Estimation
library(MCMCpack) #Needed for Long AR Regression for Estimating Innovations
library(bayesreg) #Needed for Bayesian Horseshoe/Bayesian Horseshoe+ Estimation
library(forecast) #Needed for Producing Forecasts
library(datasets) #Needed to Load CO2 Data for Mauna Loa, HI, US

options(scipen=999)
###############################################################################


###################################################
#Function Required for Obtaining Lagged Time Series
#
#Arguments: x = time series
#           k = lag
###################################################
lag.func<-function(x,k=1){
  t=length(x)
  y=c(rep(NA,t))
  for(i in (k+1):t){
    y[i]=x[i-k]
  }
  return(y)
}
###################################################



###########################################################################
#Data Used In Illustrations (Monthly CO2 Measurements in Mauna Loa, Hawaii)
###########################################################################

#Obtain Dataset and Apply Seasonal and Single Lag Differencing
maunaloa.co2=as.numeric(co2)
maunaloa.co2.time=as.numeric(time(co2))
maunaloa.co2.seasdiff.co2=c(rep(NA,12),diff(maunaloa.co2,12))
maunaloa.co2.final=c(rep(NA,1),diff(maunaloa.co2.seasdiff.co2,1))

#Split Into Modeling and Validation Periods (Mauna Loa)
MODEL.PERIOD=maunaloa.co2.time<1990
VALIDATION.PERIOD=maunaloa.co2.time>=1990

#Plot of Final Stationary Time Series Used in Illustration 
plot(x=maunaloa.co2.time,y=maunaloa.co2.final,xlab="",ylab="",type="n",
     main="CO2 After Seasonal and Regular Differencing")
points(x=maunaloa.co2.time[MODEL.PERIOD],
       y=maunaloa.co2.final[MODEL.PERIOD],
       type="l",col="black")
points(x=maunaloa.co2.time[VALIDATION.PERIOD],
       y=maunaloa.co2.final[VALIDATION.PERIOD],
       type="l",col="black",lty=3)

#Remove NA's and Subset Data For Model Fitting Period and
#Model Validation Period

#Data Fit
maunaloa.co2.train=as.numeric(na.omit(maunaloa.co2.final[MODEL.PERIOD])) 
#Forecast
maunaloa.co2.val=as.numeric(na.omit(maunaloa.co2.final[VALIDATION.PERIOD])) 
###########################################################################



###########################################################################
#Functions to Evaluate Subset Model Selection of Various Methods
#Requires Knowing True Coefficients of Data Generating Process
#Used in Simulation Experiments For Evaluation of Subset ARMA Selection
#
#Arguments: truecoef = True Known Coefficients
#           estcoef  = Estimated Coefficients
###########################################################################

#Identifies if Non-Zero Parameters Have Been Selected in Final Model (C)
id.sig.coef.func<-function(truecoef,estcoef){
  all(which(truecoef!=0) %in% which(estcoef!=0))
}

#Identifies if the True Model Has Been Selected (I)
id.true.coef.func<-function(truecoef,estcoef){
  if(length(which(truecoef!=0))==length(which(estcoef!=0))){
    all(sort(which(truecoef!=0))==as.numeric(sort(which(estcoef!=0))))
  }else{
    FALSE
  }
}

#Identifies the Proportion of Truly Non-zero Parameters Missed in Final Model (-)
false.neg.func<-function(truecoef,estcoef){
  mean(estcoef[which(truecoef!=0)]==0)
}

#Identifies the Proportion of Truly Zero Parameters Selected in Final Model (+)
false.pos.func<-function(truecoef,estcoef){
  mean(estcoef[which(truecoef==0)]!=0)
}

#Outputs a Vector of Summary Statistics (C,I,-,+)
fulleval.func<-function(truecoef,estcoef){
  return(c(id.sig.coef.func(truecoef,estcoef), #C
           id.true.coef.func(truecoef,estcoef), #I
           false.neg.func(truecoef,estcoef), #-
           false.pos.func(truecoef,estcoef))) #+
}
###########################################################################



###########################################################################
#Function to Perform ADLASSO or ADENET Subse ARMA Estimation With
#Tuning Parameter Selection Based on Minimization of AIC or BIC
#
#Arguments: x = Time Series to Be Modeled Using subset ARMA(maxP,maxQ)
#           h = Horizon Specific Model (Defaults to 1)
#           long.ar.select = Indicator Determining if Model Selection
#                            Should Be Performed in the Initial Modeling 
#                            of Long AR Process (Defaults to F)
#           maxP = Maximum AR Order 
#           maxQ = Maximum MA Order
#           updateMA = Indicator Determining if Moving Average Terms Should  
#                      Be Updated After Initial Coefficients Selected
#                      (Defaults to F)
#           BIC1 = Indicator Determining if BIC Should Be Used in Stage 1 
#                  Estimation of Weights (Defaults to T)
#           BIC2 = Indicator Determining if BIC Should Be Used in Stage 2
#                  Final Model Selection (Defaults to T)
#           alpha = Elastic Net Mixing Parameter 
#                   (0=Ridge,1=Lasso,Other=Elastic Net)
#                   (Defaults to Sequence 0,0.1,0.2,...,1)
#           eta = Exponent Applied to Weights (Defaults to 2)
#           Method = Choose Either "ADLASSO" or "ADENET"
#
#Source: Wang and Leng(2007) and Efron et al.(2004) and Zou and Hastie(2005)
###########################################################################

#Creation of Function
adshrink123.func<-function(x,h=1,long.ar.select=F,maxP,maxQ,updateMA=F,
                           BIC1=T,BIC2=T,eta=2,alpha=seq(0,1,0.1),
                           Method=c("ADLASSO","ADENET")){
  
  #Package Required
  require(glmnet) #Performs Ridge, Lasso, Elastic Net Estimation
  
  Method=match.arg(Method)
  
  Nt=length(x) #Length of Input Time Series
  
  #Fit Long AR Model to Estimate Innovations
  max.ar.order=ceiling(10*log10(Nt)) #Maximum Autoregressive Order
  init.mod.est=ar(x,aic=long.ar.select, #Allows Stepwise Selection
               order.max=max.ar.order,demean=T)
  init.mod.error=residuals(init.mod.est)
  init.mod.order=length(which(is.na(init.mod.error)))
  
  #Create Model Matrix of AR and MA terms
  dataP=foreach(p=1:maxP,.combine=cbind)%do%{
    lag.func(x,k=(p+h-1))
  }
  dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
    lag.func(init.mod.error,k=(q+h-1))
  }
  first.modX=as.matrix(cbind(dataP,
                dataQ))[-(1:(init.mod.order+max(maxP,maxQ)+h-1)),]
  first.y=x[-(1:(init.mod.order+max(maxP,maxQ)+h-1))]
  
  #Number of Alphas For Elastic Net
  n.alpha=length(alpha)
  
  #Estimation Via ADLASSO
  if(Method=="ADLASSO"){
    #Initial LASSO Weights
    first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,alpha=1)
    first.mod.RSS=colSums((first.y-predict(first.mod.est,
                                           newx=first.modX))^2)
    if(BIC1){
      first.bic.out=log(length(first.y))*first.mod.est$df+
        length(first.y)*log(as.vector(first.mod.RSS)/length(first.y))
      first.mod.lambda=first.mod.est$lambda[which.min(first.bic.out)]
    }else{
      first.aic.out=2*first.mod.est$df+
        length(first.y)*log(as.vector(first.mod.RSS)/length(first.y))
      first.mod.lambda=first.mod.est$lambda[which.min(first.aic.out)]
    }
    first.mod.coef=as.numeric(coef(first.mod.est,
                   s=first.mod.lambda,method="lambda"))[-1]
    first.mod.mu=as.numeric(coef(first.mod.est,
                                 s=first.mod.lambda,method="lambda"))[1]
    weights=abs(first.mod.coef+1/length(first.y))^(-eta)
    
    #Update Model Matrix of AR and MA terms Based off Initial Estimation
    if(updateMA){
      update.mod.predict=rep(NA,length(x))
      update.mod.error=rep(0,length(x))
      for(v in (h+max(maxP,maxQ)):Nt){
        update.mod.predict[v]=first.mod.mu+
          x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
          update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
        update.mod.error[v]=x[v]-update.mod.predict[v]
      }
      update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
        lag.func(update.mod.error,k=(q+h-1))
      }
      second.modX=as.matrix(cbind(dataP,
                      update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
      second.y=x[-(1:(max(maxP,maxQ)+h-1))]
    }else{
      second.modX=first.modX
      second.y=first.y
    }
    
    second.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,alpha=1,
                          thresh=1e-16, penalty.factor=weights)
    second.mod.RSS=colSums((second.y-predict(second.mod.est,
                                             newx=second.modX))^2)
    if(BIC2){
      second.bic.out=log(length(second.y))*second.mod.est$df+
        length(second.y)*log(as.vector(second.mod.RSS)/length(second.y))
      second.mod.lambda=second.mod.est$lambda[which.min(second.bic.out)]
    }else{
      second.aic.out=2*second.mod.est$df+
        length(second.y)*log(as.vector(second.mod.RSS)/length(second.y))
      second.mod.lambda=second.mod.est$lambda[which.min(second.aic.out)]
    }
    
    final.mod.coef=as.numeric(coef(second.mod.est,s=second.mod.lambda,
                   method="lambda"))[-1]
    nonzero.select=which(final.mod.coef!=0)
    final.mod.int=as.numeric(coef(second.mod.est,s=second.mod.lambda,
                  method="lambda"))[1]
    final.mod.s2=sum((second.y-predict(second.mod.est,newx=second.modX,
                 s=second.mod.lambda,method="lambda"))^2)/(length(second.y)-
                 sum(final.mod.coef[nonzero.select]!=0)-1)
    
    out=list(final.mod.coef=final.mod.coef, #Final Selection of Coefficients
             final.mod.int=final.mod.int,   #Final Estimated Intercept
             final.mod.s2=final.mod.s2,     #Final Estimated Noise Variance
             nonzero.select=nonzero.select) #Identifies the Nonzero Parameters
  }
  
  #Estimation Via ADENET
  if(Method=="ADENET"){
    #Initial LASSO Weights
    first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,alpha=1)
    first.mod.RSS=colSums((first.y-predict(first.mod.est,newx=first.modX))^2)
    if(BIC1){
      first.bic.out=log(length(first.y))*first.mod.est$df+
        length(first.y)*log(as.vector(first.mod.RSS)/length(first.y))
      first.mod.lambda=first.mod.est$lambda[which.min(first.bic.out)]
      first.cv.out=c(1,first.mod.lambda,min(first.bic.out))
    }else{
      first.aic.out=2*first.mod.est$df+
        length(first.y)*log(as.vector(first.mod.RSS)/length(first.y))
      first.mod.lambda=first.mod.est$lambda[which.min(first.aic.out)]
      first.cv.out=c(1,first.mod.lambda,min(first.aic.out))
    }
    
    first.mod.alpha=1
    first.mod.lambda=first.cv.out[2]
    first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,
                  alpha=first.mod.alpha,lambda=first.mod.lambda)
    first.mod.coef=as.numeric(coef(first.mod.est))[-1]
    first.mod.mu=as.numeric(coef(first.mod.est))[1]
    
    weights=abs(first.mod.coef+1/length(first.y))^(-eta)
    
    #Update Model Matrix of AR and MA terms Based off Initial Estimation
    if(updateMA){
      update.mod.predict=rep(NA,length(x))
      update.mod.error=rep(0,length(x))
      for(v in (h+max(maxP,maxQ)):Nt){
        update.mod.predict[v]=first.mod.mu+
          x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
          update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
        update.mod.error[v]=x[v]-update.mod.predict[v]
      }
      update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
        lag.func(update.mod.error,k=(q+h-1))
      }
      second.modX=as.matrix(cbind(dataP,
                                  update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
      second.y=x[-(1:(max(maxP,maxQ)+h-1))]
    }else{
      second.modX=first.modX
      second.y=first.y
    }
    
    #Final Elastic Net Estimates (Search Through All Lambdas)
    second.cv.out=foreach(a=1:n.alpha,.combine=rbind)%do%{
      second.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,
                     alpha=alpha[a],penalty.factor=weights)
      second.mod.RSS=colSums((second.y-predict(second.mod.est,
                                               newx=second.modX))^2)
      if(BIC2){
        second.bic.out=log(length(second.y))*second.mod.est$df+
          length(second.y)*log(as.vector(second.mod.RSS)/length(second.y))
        second.mod.lambda=second.mod.est$lambda[which.min(second.bic.out)]
        result=c(alpha[a],second.mod.lambda,min(second.bic.out))
      }else{
        second.aic.out=2*second.mod.est$df+
          length(second.y)*log(as.vector(second.mod.RSS)/length(second.y))
        second.mod.lambda=second.mod.est$lambda[which.min(second.aic.out)]
        result=c(alpha[a],second.mod.lambda,min(second.aic.out))
      }
      result
    }  
    
    second.mod.alpha=alpha[which.min(second.cv.out[,3])]
    second.mod.lambda=second.cv.out[which.min(second.cv.out[,3]),2]
    second.mod.est=glmnet(y=second.y,x=second.modX,
                          standardize=T,alpha=second.mod.alpha,
                          lambda=second.mod.lambda,penalty.factor=weights)
    second.mod.coef=as.numeric(coef(second.mod.est))[-1]
    second.mod.mu=as.numeric(coef(second.mod.est))[1]
    
    final.mod.coef=second.mod.coef
    nonzero.select=which(final.mod.coef!=0)
    final.mod.int=second.mod.mu
    final.mod.s2=sum((second.y-predict(second.mod.est,
                      newx=second.modX))^2)/(length(second.y)-
                      sum(final.mod.coef[nonzero.select]!=0)-1)
    
    out=list(final.mod.coef=final.mod.coef, #Final Selection of Coefficients
             final.mod.int=final.mod.int,   #Final Estimated Intercept
             final.mod.s2=final.mod.s2,     #Final Estimated Noise Variance
             nonzero.select=nonzero.select) #Identifies the Nonzero Parameters
  }
  
  return(out)
}

#Illustration of Function for ADLASSO Estimation
adlasso1=adshrink123.func(x=maunaloa.co2.train,h=1,
         long.ar.select=F,maxP=14,maxQ=14,
         updateMA=F,BIC1=F,BIC2=F,eta=2,
         alpha=seq(0,1,0.1),Method="ADLASSO")
adlasso2=adshrink123.func(x=maunaloa.co2.train,h=1,
         long.ar.select=F,maxP=14,maxQ=14,
         updateMA=F,BIC1=F,BIC2=T,eta=2,
         alpha=seq(0,1,0.1),Method="ADLASSO")
adlasso3=adshrink123.func(x=maunaloa.co2.train,h=1,
         long.ar.select=F,maxP=14,maxQ=14,
         updateMA=F,BIC1=T,BIC2=T,eta=2,
         alpha=seq(0,1,0.1),Method="ADLASSO")

#Illustration of Function for ADENET Estimation
adenet1=adshrink123.func(x=maunaloa.co2.train,h=1,
        long.ar.select=F,maxP=14,maxQ=14,
        updateMA=F,BIC1=F,BIC2=F,eta=2,
        alpha=seq(0,1,0.1),Method="ADENET")
adenet2=adshrink123.func(x=maunaloa.co2.train,h=1,
        long.ar.select=F,maxP=14,maxQ=14,
        updateMA=F,BIC1=F,BIC2=T,eta=2,
        alpha=seq(0,1,0.1),Method="ADENET")
adenet3=adshrink123.func(x=maunaloa.co2.train,h=1,
        long.ar.select=F,maxP=14,maxQ=14,
        updateMA=F,BIC1=T,BIC2=T,eta=2,
        alpha=seq(0,1,0.1),Method="ADENET")
###########################################################################



###########################################################################
#Functions Used  to Partition Data into a Train Set containing the 
#   first (1-test.per)x100% of data and a Test Set containing the 
#   last (test.per)x100% of data. The second function 
#   institutes a gap of length max.pq to remove temporal dependencies 
#Arguments: x = time series
#           max.pq = maximum ar/ma order to consider 
#                   (covers temporal dependence)
#           test.per = Percent of Data to Consider in Test Set
#                      (Defaults to 0.2)
#
#Output: Vector of O's (Fit) and 1's (Tuning Parameter Selection)
###########################################################################

#Function Splitting Data According to Classic Out-of-Sample Procedure
OOS.IndepCV.func<-function(x,test.per=0.20){
  N=length(x)
  test.N=ceiling(test.per*N)
  if(test.per>0.5){ 
    warning("Less than Half the Dataset Is Being Used for Training ")
  }
  Block.Vector=rep(0,N)
  Block.Vector[(N-test.N+1):N]=1
  return(as.matrix(Block.Vector))
}

#Function Splitting OOS But Removing Last max.pq from End of Training Set
OOS.DepCV.func<-function(x,max.pq=NULL,test.per=0.20){
  N=length(x)
  test.N=ceiling(test.per*N)
  if(is.null(max.pq)) max.pq=ceiling(10*log10(test.N))
  if(test.per>0.5){
    warning("Less than Half the Dataset Is Being Used for Training ")
  }
  if(max.pq>(N-test.N)) warning("max.pq > Number of Obs in Training")
  Block.Vector=rep(0,N)
  Block.Vector[(N-test.N+1):N]=1
  Block.Vector[(N-test.N-max.pq+1):(N-test.N)]=NA
  return(as.matrix(Block.Vector))
}
###########################################################################



###########################################################################
#Function to Perform ADLASSO or ADENET Subse ARMA Estimation With
#Tuning Parameter Selection Based on Out-Of-Sample Period
#
#Arguments: x = Time Series to Be Modeled Using subset ARMA(maxP,maxQ)
#           h = Horizon Specific Model (Defaults to 1)
#           long.ar.select = Indicator Determining if Model Selection Is 
#                            Performed in the Initial Modeling of  
#                            Long AR Process (Defaults to F)
#           maxP = Maximum AR Order 
#           maxQ = Maximum MA Order
#           test.per = Percentage of Data Removed for 
#                      Tuning Parameter Selection (Defaults to 0.2 => 20%)
#           max.pq = Maximimum Temporal Dependence Considered in depOOS
#                    (Defaults to max(maxP,maxQ))
#           updateMA = Indicator Determining if Moving Average Terms Is 
#                      Updated After Initial Coefficients Selected
#                      (Defaults to F)
#           BIC1 = Indicator Determining if BIC Should Be Used in Stage 1 
#                  Estimation of Weights (Defaults to T)
#           BIC2 = Indicator Determining if BIC Should Be Used in Stage 2
#                  Final Model Selection (Defaults to T)
#           alpha = Elastic Net Mixing Parameter 
#                   (0=Ridge,1=Lasso,Other=Elastic Net)
#                   (Defaults to Sequence 0,0.1,0.2,...,1)
#           eta = Exponent Applied to Weights (Defaults to 2)
#           Method = Choose Either "ADLASSO" or "ADENET"
#           CV = Choose Either "OOS" or "depOOS
#           ADENET.final = Appoach for Choosing Final Lambda for Each Alpha
#                          "min" = Choose Based on Minimum MSE
#                          "1se" = Choose Largest Lambda that Results in MSE
#                                  within 1 Standard Error of the Minimum
#
#Source: Wang and Leng(2007) and Efron et al.(2004) and Zou and Hastie(2005)
###########################################################################

#Creation of Function
adshrink45.func<-function(x,h=1,long.ar.select=F,maxP,maxQ,updateMA=F,
                  BIC1=T,BIC2=T,eta=2,alpha=seq(0,1,0.1),test.per=0.2,
                  max.pq = max(maxP,maxQ),ADENET.final=c("min","1se"),
                  Method=c("ADLASSO","ADENET"),CV=c("OOS","depOOS")){
  
  #Package Required
  require(glmnet) #Performs Ridge, Lasso, Elastic Net Estimation
  
  Method=match.arg(Method)
  CV=match.arg(CV)
  ADENET.final=match.arg(ADENET.final)
  
  Nt=length(x) #Length of Input Time Series
  
  #Fit Long AR Model to Estimate Innovations
  max.ar.order=ceiling(10*log10(Nt)) #Maximum Autoregressive Order
  init.mod.est=ar(x,aic=long.ar.select, #Allows Stepwise Selection
                  order.max=max.ar.order,demean=T)
  init.mod.error=residuals(init.mod.est)
  init.mod.order=length(which(is.na(init.mod.error)))
  
  #Create Model Matrix of AR and MA terms
  dataP=foreach(p=1:maxP,.combine=cbind)%do%{
    lag.func(x,k=(p+h-1))
  }
  dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
    lag.func(init.mod.error,k=(q+h-1))
  }
  first.modX=as.matrix(cbind(dataP,
                  dataQ))[-(1:(init.mod.order+max(maxP,maxQ)+h-1)),]
  first.y=x[-(1:(init.mod.order+max(maxP,maxQ)+h-1))]
  
  #Number of Alphas For Elastic Net
  n.alpha=length(alpha)
  
  #Estimation Via ADLASSO
  if(Method=="ADLASSO"){
    if(CV=="OOS"){
      Fold.Vector=OOS.IndepCV.func(x=first.y,test.per=test.per)
      
      in.train=which(Fold.Vector==0)
      in.test=which(Fold.Vector==1)
      
      first.mod.est=glmnet(y=first.y[in.train],x=first.modX[in.train,],
                    standardize=T,alpha=1)
      first.mod.res=(first.y[in.test]-predict(first.mod.est,
                    newx=first.modX[in.test,]))^2
      
      OOS.MSE1=apply(first.mod.res,2,mean)
      lambda1.min=first.mod.est$lambda[which.min(OOS.MSE1)]
      lambda1.1se=first.mod.est$lambda[min(which(OOS.MSE1<(min(OOS.MSE1)+
                  sd(OOS.MSE1)/sqrt(length(OOS.MSE1)))))]
      
      first.mod.coef=as.numeric(coef(first.mod.est,s=lambda1.min))[-1]
      first.mod.mu=as.numeric(coef(first.mod.est,s=lambda1.min))[1]
      weights=abs(first.mod.coef+1/length(first.y))^(-eta)
      
      if(updateMA){
        update.mod.predict=rep(NA,length(x))
        update.mod.error=rep(0,length(x))
        for(v in (h+max(maxP,maxQ)):Nt){
          update.mod.predict[v]=first.mod.mu+
            x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
            update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
          update.mod.error[v]=x[v]-update.mod.predict[v]
        }
        update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
          lag.func(update.mod.error,k=(q+h-1))
        }
        second.modX=as.matrix(cbind(dataP,
                          update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
        second.y=x[-(1:(max(maxP,maxQ)+h-1))]
      }else{
        second.modX=first.modX
        second.y=first.y
      }
      
      second.mod.est=glmnet(y=second.y[in.train],x=second.modX[in.train,],
                     standardize=T,alpha=1,penalty.factor=weights)
      second.mod.res=(second.y[in.test]-
                     predict(second.mod.est,newx=second.modX[in.test,]))^2
      
      OOS.MSE2=apply(second.mod.res,2,mean)
      lambda2.min=second.mod.est$lambda[which.min(OOS.MSE2)]
      lambda2.1se=second.mod.est$lambda[min(which(OOS.MSE2<(min(OOS.MSE2)+
                  sd(OOS.MSE2)/sqrt(length(OOS.MSE2)))))]
      
      second.mod.coef=as.numeric(coef(second.mod.est,s=lambda2.1se))[-1]
      second.mod.mu=as.numeric(coef(second.mod.est,s=lambda2.1se))[1]
      
      final.mod.coef=second.mod.coef
      nonzero.select=which(final.mod.coef!=0)
      final.mod.int=second.mod.mu
      final.mod.s2=sum((second.y-predict(second.mod.est,newx=second.modX,
                   s=lambda2.1se,method="lambda"))^2)/(length(second.y)-
                  sum(final.mod.coef[nonzero.select]!=0)-1)
      
      out=list(final.mod.coef=final.mod.coef, #Final Selected of Coefficients
               final.mod.int=final.mod.int,   #Final Estimated Intercept
               final.mod.s2=final.mod.s2,     #Final Estimated Noise Variance
               nonzero.select=nonzero.select) #Identifies Nonzero Parameters 
    }
    
    if(CV=="depOOS"){
      Fold.Vector=OOS.DepCV.func(x=first.y,test.per=test.per)
      
      in.train=which(Fold.Vector==0)
      in.test=which(Fold.Vector==1)
      
      first.mod.est=glmnet(y=first.y[in.train],
                    x=first.modX[in.train,],standardize=T,alpha=1)
      first.mod.res=(first.y[in.test]-predict(first.mod.est,
                    newx=first.modX[in.test,]))^2
      
      OOS.MSE1=apply(first.mod.res,2,mean)
      lambda1.min=first.mod.est$lambda[which.min(OOS.MSE1)]
      lambda1.1se=first.mod.est$lambda[min(which(OOS.MSE1<(min(OOS.MSE1)+
                  sd(OOS.MSE1)/sqrt(length(OOS.MSE1)))))]
      
      first.mod.coef=as.numeric(coef(first.mod.est,s=lambda1.min))[-1]
      first.mod.mu=as.numeric(coef(first.mod.est,s=lambda1.min))[1]
      weights=abs(first.mod.coef+1/length(first.y))^(-eta)
      
      if(updateMA){
        update.mod.predict=rep(NA,length(x))
        update.mod.error=rep(0,length(x))
        for(v in (h+max(maxP,maxQ)):Nt){
          update.mod.predict[v]=first.mod.mu+
            x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
            update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
          update.mod.error[v]=x[v]-update.mod.predict[v]
        }
        update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
          lag.func(update.mod.error,k=(q+h-1))
        }
        second.modX=as.matrix(cbind(dataP,
                          update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
        second.y=x[-(1:(max(maxP,maxQ)+h-1))]
      }else{
        second.modX=first.modX
        second.y=first.y
      }
      
      second.mod.est=glmnet(y=second.y[in.train],x=second.modX[in.train,],
                     standardize=T,alpha=1,penalty.factor=weights)
      second.mod.res=(second.y[in.test]-predict(second.mod.est,
                     newx=second.modX[in.test,]))^2
      
      OOS.MSE2=apply(second.mod.res,2,mean)
      lambda2.min=second.mod.est$lambda[which.min(OOS.MSE2)]
      lambda2.1se=second.mod.est$lambda[min(which(OOS.MSE2<(min(OOS.MSE2)+
                  sd(OOS.MSE2)/sqrt(length(OOS.MSE2)))))]
      
      second.mod.coef=as.numeric(coef(second.mod.est,s=lambda2.1se))[-1]
      second.mod.mu=as.numeric(coef(second.mod.est,s=lambda2.1se))[1]
      
      final.mod.coef=second.mod.coef
      nonzero.select=which(final.mod.coef!=0)
      final.mod.int=second.mod.mu
      final.mod.s2=sum((second.y-predict(second.mod.est,newx=second.modX,
                   s=lambda2.1se,method="lambda"))^2)/
                   (length(second.y)-sum(final.mod.coef[nonzero.select]!=0)-1)
      
      out=list(final.mod.coef=final.mod.coef, #Final Selected of Coefficients
               final.mod.int=final.mod.int,   #Final Estimated Intercept
               final.mod.s2=final.mod.s2,     #Final Estimated Noise Variance
               nonzero.select=nonzero.select) #Identifies Nonzero Parameters   
      } 
  }
  
  #Estimation Via ADENET
  if(Method=="ADENET"){
    if(CV=="OOS"){
      Fold.Vector=OOS.IndepCV.func(x=first.y,test.per=test.per)
      
      in.train=which(Fold.Vector==0)
      in.test=which(Fold.Vector==1)
      
      first.mod.est=glmnet(y=first.y[in.train],
                           x=first.modX[in.train,],
                           standardize=T,alpha=1)
      first.mod.res=(first.y[in.test]-
                    predict(first.mod.est,newx=first.modX[in.test,]))^2
      OOS.MSE1=apply(first.mod.res,2,mean)
      lambda1.min=first.mod.est$lambda[which.min(OOS.MSE1)]
      lambda1.1se=first.mod.est$lambda[min(which(OOS.MSE1<(min(OOS.MSE1)+
                  sd(OOS.MSE1)/sqrt(length(OOS.MSE1)))))]
      first.cv.out=c(1,lambda1.min,OOS.MSE1[which.min(OOS.MSE1)])
      
      first.mod.alpha=1
      first.mod.lambda=first.cv.out[2]
      first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,
                    alpha=first.mod.alpha,lambda=first.mod.lambda)
      first.mod.coef=as.numeric(coef(first.mod.est))[-1]
      first.mod.mu=as.numeric(coef(first.mod.est))[1]
      
      weights=abs(first.mod.coef+1/length(first.y))^(-eta)
      
      if(updateMA){
        update.mod.predict=rep(NA,length(x))
        update.mod.error=rep(0,length(x))
        for(v in (h+max(maxP,maxQ)):Nt){
          update.mod.predict[v]=first.mod.mu+
            x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
            update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
          update.mod.error[v]=x[v]-update.mod.predict[v]
        }
        update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
          lag.func(update.mod.error,k=(q+h-1))
        }
        second.modX=as.matrix(cbind(dataP,
                          update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
        second.y=x[-(1:(max(maxP,maxQ)+h-1))]
      }else{
        second.modX=first.modX
        second.y=first.y
      }
      
      second.cv.out=foreach(a=1:n.alpha,.combine=rbind)%do%{
        second.mod.est=glmnet(y=second.y[in.train],
                              x=second.modX[in.train,],standardize=T,
                              alpha=alpha[a])
        second.mod.res=(second.y[in.test]-
                       predict(second.mod.est,newx=second.modX[in.test,]))^2
        OOS.MSE2=apply(second.mod.res,2,mean)
        lambda2.min=second.mod.est$lambda[which.min(OOS.MSE2)]
        lambda2.1se=second.mod.est$lambda[min(which(OOS.MSE2<(min(OOS.MSE2)+
                    sd(OOS.MSE2)/sqrt(length(OOS.MSE2)))))]
        if(ADENET.final=="min") out=c(alpha[a],lambda2.min,min(OOS.MSE2))
        if(ADENET.final=="1se") out=c(alpha[a],lambda2.1se,
                OOS.MSE2[min(which(OOS.MSE2<(min(OOS.MSE2)+
                sd(OOS.MSE2)/sqrt(length(OOS.MSE2)))))])
        out
      }
      
      second.mod.alpha=alpha[which.min(second.cv.out[,3])]
      second.mod.lambda=second.cv.out[which.min(second.cv.out[,3]),2]
      second.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,
                     alpha=second.mod.alpha,lambda=second.mod.lambda,
                     penalty.factor=weights)
      second.mod.coef=as.numeric(coef(second.mod.est))[-1]
      second.mod.mu=as.numeric(coef(second.mod.est))[1]
      
      final.mod.coef=second.mod.coef
      nonzero.select=which(final.mod.coef!=0)
      final.mod.int=second.mod.mu
      final.mod.s2=sum((second.y-
                   predict(second.mod.est,
                   newx=second.modX))^2)/(length(second.y)-
                   sum(final.mod.coef[nonzero.select]!=0)-1)
      
      out=list(final.mod.coef=final.mod.coef, #Final Selected of Coefficients
               final.mod.int=final.mod.int,   #Final Estimated Intercept
               final.mod.s2=final.mod.s2,     #Final Estimated Noise Variance
               nonzero.select=nonzero.select) #Identifies Nonzero Parameters 
    }
  
    if(CV=="depOOS"){
      Fold.Vector=OOS.DepCV.func(x=first.y,test.per=test.per)
      
      in.train=which(Fold.Vector==0)
      in.test=which(Fold.Vector==1)
      
      first.mod.est=glmnet(y=first.y[in.train],
                           x=first.modX[in.train,],
                           standardize=T,alpha=1)
      first.mod.res=(first.y[in.test]-
                    predict(first.mod.est,newx=first.modX[in.test,]))^2
      OOS.MSE1=apply(first.mod.res,2,mean)
      lambda1.min=first.mod.est$lambda[which.min(OOS.MSE1)]
      lambda1.1se=first.mod.est$lambda[min(which(OOS.MSE1<(min(OOS.MSE1)+
                  sd(OOS.MSE1)/sqrt(length(OOS.MSE1)))))]
      first.cv.out=c(1,lambda1.min,OOS.MSE1[which.min(OOS.MSE1)])
      
      first.mod.alpha=1
      first.mod.lambda=first.cv.out[2]
      first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,
                    alpha=first.mod.alpha,lambda=first.mod.lambda)
      first.mod.coef=as.numeric(coef(first.mod.est))[-1]
      first.mod.mu=as.numeric(coef(first.mod.est))[1]
      
      weights=abs(first.mod.coef+1/length(first.y))^(-eta)
      
      if(updateMA){
        update.mod.predict=rep(NA,length(x))
        update.mod.error=rep(0,length(x))
        for(v in (h+max(maxP,maxQ)):Nt){
          update.mod.predict[v]=first.mod.mu+
            x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
            update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
          update.mod.error[v]=x[v]-update.mod.predict[v]
        }
        update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
          lag.func(update.mod.error,k=(q+h-1))
        }
        second.modX=as.matrix(cbind(dataP,
                          update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
        second.y=x[-(1:(max(maxP,maxQ)+h-1))]
      }else{
        second.modX=first.modX
        second.y=first.y
      }
      
      second.cv.out=foreach(a=1:n.alpha,.combine=rbind)%do%{
        second.mod.est=glmnet(y=second.y[in.train],
                              x=second.modX[in.train,],standardize=T,
                              alpha=alpha[a])
        second.mod.res=(second.y[in.test]-
                       predict(second.mod.est,
                       newx=second.modX[in.test,]))^2
        OOS.MSE2=apply(second.mod.res,2,mean)
        lambda2.min=second.mod.est$lambda[which.min(OOS.MSE2)]
        lambda2.1se=second.mod.est$lambda[min(which(OOS.MSE2<(min(OOS.MSE2)+
                    sd(OOS.MSE2)/sqrt(length(OOS.MSE2)))))]
        if(ADENET.final=="min") out=c(alpha[a],lambda2.min,min(OOS.MSE2))
        if(ADENET.final=="1se") out=c(alpha[a],lambda2.1se,
                  OOS.MSE2[min(which(OOS.MSE2<(min(OOS.MSE2)+
                  sd(OOS.MSE2)/sqrt(length(OOS.MSE2)))))])
        out
      }
      
      second.mod.alpha=alpha[which.min(second.cv.out[,3])]
      second.mod.lambda=second.cv.out[which.min(second.cv.out[,3]),2]
      second.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,
                     alpha=second.mod.alpha,lambda=second.mod.lambda,
                     penalty.factor=weights)
      second.mod.coef=as.numeric(coef(second.mod.est))[-1]
      second.mod.mu=as.numeric(coef(second.mod.est))[1]
      
      final.mod.coef=second.mod.coef
      nonzero.select=which(final.mod.coef!=0)
      final.mod.int=second.mod.mu
      final.mod.s2=sum((second.y-
                   predict(second.mod.est,
                   newx=second.modX))^2)/(length(second.y)-
                   sum(final.mod.coef[nonzero.select]!=0)-1)
      
      out=list(final.mod.coef=final.mod.coef, #Final Selected of Coefficients
               final.mod.int=final.mod.int,   #Final Estimated Intercept
               final.mod.s2=final.mod.s2,     #Final Estimated Noise Variance
               nonzero.select=nonzero.select) #Identifies Nonzero Parameters 
    }
  }
  
  return(out)
}

#Illustration of Function for ADLASSO Estimation
adlasso4=adshrink45.func(x=maunaloa.co2.train,h=1,
         long.ar.select=F,maxP=14,maxQ=14,
         updateMA=F,BIC1=F,BIC2=F,eta=2,alpha=seq(0,1,0.1),Method="ADLASSO",
         test.per=0.2,CV="OOS")
adlasso5=adshrink45.func(x=maunaloa.co2.train,h=1,
         long.ar.select=F,maxP=14,maxQ=14,
         updateMA=F,BIC1=F,BIC2=T,eta=2,alpha=seq(0,1,0.1),Method="ADLASSO",
         test.per=0.2,max.pq=max(maxP,maxQ),CV="depOOS")

#Illustration of Function for ADENET Estimation
adenet4=adshrink45.func(x=maunaloa.co2.train,h=1,
        long.ar.select=F,maxP=14,maxQ=14,
        updateMA=F,BIC1=F,BIC2=F,eta=2,alpha=seq(0,1,0.1),Method="ADENET",
        ADENET.final="min",test.per=0.2,CV="OOS")
adenet5=adshrink45.func(x=maunaloa.co2.train,h=1,
        long.ar.select=F,maxP=14,maxQ=14,
        updateMA=F,BIC1=F,BIC2=T,eta=2,alpha=seq(0,1,0.1),Method="ADENET",
        ADENET.final="min",test.per=0.2,max.pq=max(maxP,maxQ),CV="depOOS")
##############################################################################



##############################################################################
#Function to Perform ADLASSO or ADENET Subse ARMA Estimation With
#Tuning Parameter Selection Based on Regular K-fold Cross Validation 
#
#Arguments: x = Time Series to Be Modeled Using subset ARMA(maxP,maxQ)
#           h = Horizon Specific Model (Defaults to 1)
#           long.ar.select = Indicator Determining if Model Selection Is 
#                            Performed in the Initial Modeling of  
#                            Long AR Process (Defaults to F)
#           maxP = Maximum AR Order 
#           maxQ = Maximum MA Order
#           updateMA = Indicator Determining if Moving Average Terms Is 
#                      Updated After Initial Coefficients Selected
#                      (Defaults to F)
#           BIC1 = Indicator Determining if BIC Should Be Used in Stage 1 
#                  Estimation of Weights (Defaults to T)
#           BIC2 = Indicator Determining if BIC Should Be Used in Stage 2
#                  Final Model Selection (Defaults to T)
#           alpha = Elastic Net Mixing Parameter 
#                   (0=Ridge,1=Lasso,Other=Elastic Net)
#                   (Defaults to Sequence 0,0.1,0.2,...,1)
#           eta = Exponent Applied to Weights (Defaults to 2)
#           Method = Choose Either "ADLASSO" or "ADENET"
#           K = Number of Folds for General CV (Defaults to NULL => LOOCV)
#           ADENET.final = Appoach for Choosing Final Lambda for Each Alpha
#                          "min" = Choose Based on Minimum MSE
#                          "1se" = Choose Largest Lambda that Results in MSE
#                                  within 1 Standard Error of the Minimum
#
#Source: Wang and Leng(2007) and Efron et al.(2004) and Zou and Hastie(2005)
##############################################################################

#Creation of Function
adshrink678.func<-function(x,h=1,long.ar.select=F,maxP,maxQ,updateMA=F,
                          BIC1=T,BIC2=T,eta=2,alpha=seq(0,1,0.1),
                          ADENET.final=c("min","1se"),
                          K=NULL,Method=c("ADLASSO","ADENET")){
  
  #Package Required
  require(glmnet) #Performs Ridge, Lasso, Elastic Net Estimation
  
  Method=match.arg(Method)
  ADENET.final=match.arg(ADENET.final)
  
  Nt=length(x) #Length of Input Time Series
  
  #Fit Long AR Model to Estimate Innovations
  max.ar.order=ceiling(10*log10(Nt)) #Maximum Autoregressive Order
  init.mod.est=ar(x,aic=long.ar.select, #Allows Stepwise Selection
                  order.max=max.ar.order,demean=T)
  init.mod.error=residuals(init.mod.est)
  init.mod.order=length(which(is.na(init.mod.error)))
  
  #Create Model Matrix of AR and MA terms
  dataP=foreach(p=1:maxP,.combine=cbind)%do%{
    lag.func(x,k=(p+h-1))
  }
  dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
    lag.func(init.mod.error,k=(q+h-1))
  }
  first.modX=as.matrix(cbind(dataP,
             dataQ))[-(1:(init.mod.order+max(maxP,maxQ)+h-1)),]
  first.y=x[-(1:(init.mod.order+max(maxP,maxQ)+h-1))]
  
  #Number of Alphas For Elastic Net
  n.alpha=length(alpha)
  
  #Estimation Via ADLASSO
  if(Method=="ADLASSO"){
    if(is.null(K)){
      first.mod.est=cv.glmnet(y=first.y,x=first.modX,standardize=T,
                    alpha=1,parallel=F)
    }else{
      first.mod.est=cv.glmnet(y=first.y,x=first.modX,standardize=T,
                    alpha=1,nfolds=K,parallel=F)
    }
    
    first.mod.coef=coef(first.mod.est,s=first.mod.est$lambda.min)[-1]
    first.mod.mu=coef(first.mod.est,s=first.mod.est$lambda.min)[1]
    weights=abs(first.mod.coef+1/length(first.y))^(-2)
    
    if(updateMA){
      update.mod.predict=rep(NA,length(x))
      update.mod.error=rep(0,length(x))
      for(v in (h+max(maxP,maxQ)):Nt){
        update.mod.predict[v]=first.mod.mu+
          x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
          update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
        update.mod.error[v]=x[v]-update.mod.predict[v]
      }
      update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
        lag.func(update.mod.error,k=(q+h-1))
      }
      second.modX=as.matrix(cbind(dataP,
                  update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
      second.y=x[-(1:(max(maxP,maxQ)+h-1))]
    }else{
      second.modX=first.modX
      second.y=first.y
    }
    
    if(is.null(K)){
      second.mod.est=cv.glmnet(y=second.y,x=second.modX,standardize=T,
                     parallel=F,alpha=1,penalty.factor=weights)
    }else{
      second.mod.est=cv.glmnet(y=second.y,x=second.modX,standardize=T,
                     parallel=F,alpha=1,nfolds=K,penalty.factor=weights)
    }
    
    final.mod.coef=coef(second.mod.est,s=second.mod.est$lambda.1se)[-1]
    nonzero.select=which(final.mod.coef!=0)
    final.mod.int=coef(second.mod.est,s=second.mod.est$lambda.1se)[1]
    final.mod.s2=sum((second.y-predict(second.mod.est,newx=second.modX,
                 s=second.mod.est$lambda.1se,method="lambda"))^2)/
                 (length(second.y)-sum(final.mod.coef[nonzero.select]!=0)-1)
    
    out=list(final.mod.coef=final.mod.coef, #Final Selection of Coefficients
             final.mod.int=final.mod.int,   #Final Estimated Intercept
             final.mod.s2=final.mod.s2,     #Final Estimated Noise Variance
             nonzero.select=nonzero.select) #Identifies the Nonzero Parameters
  }
  
  #Estimation Via ADENET
  if(Method=="ADENET"){
    if(is.null(K)){
      first.mod.est=cv.glmnet(parallel=F,y=first.y,x=first.modX,
                    standardize=T,alpha=1)
      first.cv.out=c(1,first.mod.est$lambda.min,
        first.mod.est$cvm[which(first.mod.est$lambda==first.mod.est$lambda.min)])
    }else{
      first.mod.est=cv.glmnet(parallel=F,y=first.y,x=first.modX,
                    standardize=T,alpha=1,nfolds=K)
      first.cv.out=c(1,first.mod.est$lambda.min,
        first.mod.est$cvm[which(first.mod.est$lambda==first.mod.est$lambda.min)])
    }
    
    first.mod.alpha=1
    first.mod.lambda=first.cv.out[2]
    first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,
                  alpha=first.mod.alpha,lambda=first.mod.lambda)
    first.mod.coef=as.numeric(coef(first.mod.est))[-1]
    first.mod.mu=as.numeric(coef(first.mod.est))[1]
    
    weights=abs(first.mod.coef+1/length(first.y))^(-eta)
    
    
    if(updateMA){
      update.mod.predict=rep(NA,length(x))
      update.mod.error=rep(0,length(x))
      for(v in (h+max(maxP,maxQ)):Nt){
        update.mod.predict[v]=first.mod.mu+
          x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
          update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
        update.mod.error[v]=x[v]-update.mod.predict[v]
      }
      update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
        lag.func(update.mod.error,k=(q+h-1))
      }
      second.modX=as.matrix(cbind(dataP,
                  update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
      second.y=x[-(1:(max(maxP,maxQ)+h-1))]
    }else{
      second.modX=first.modX
      second.y=first.y
    }
    
    n.alpha=length(alpha)
    
    if(is.null(K)){
      second.cv.out=NULL
      for(a in 1:n.alpha){  
        second.mod.est=cv.glmnet(parallel=F,y=second.y,x=second.modX,
                       standardize=T,alpha=alpha[a],penalty.factor=weights)
        if(ADENET.final=="min") second.cv.out=rbind(second.cv.out,
          c(alpha[a],second.mod.est$lambda.min,
    second.mod.est$cvm[which(second.mod.est$lambda==second.mod.est$lambda.min)]))
        if(ADENET.final=="1se") second.cv.out=rbind(second.cv.out,
          c(alpha[a],second.mod.est$lambda.1se,
    second.mod.est$cvm[which(second.mod.est$lambda==second.mod.est$lambda.1se)]))
      }
    }else{
      second.cv.out=NULL
      for(a in 1:n.alpha){  
        second.mod.est=cv.glmnet(parallel=F,y=second.y,x=second.modX,
                       standardize=T,alpha=alpha[a],nfolds=K,
                       penalty.factor=weights)
        if(ADENET.final=="min") second.cv.out=rbind(second.cv.out,c(alpha[a],
            second.mod.est$lambda.min,
    second.mod.est$cvm[which(second.mod.est$lambda==second.mod.est$lambda.min)]))
        if(ADENET.final=="1se") second.cv.out=rbind(second.cv.out,
            c(alpha[a],second.mod.est$lambda.1se,
    second.mod.est$cvm[which(second.mod.est$lambda==second.mod.est$lambda.1se)]))
      }
    }
    
    second.mod.alpha=alpha[which.min(second.cv.out[,3])]
    second.mod.lambda=second.cv.out[which.min(second.cv.out[,3]),2]
    second.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,
              alpha=second.mod.alpha,lambda=second.mod.lambda,
              penalty.factor=weights)
    second.mod.coef=as.numeric(coef(second.mod.est))[-1]
    second.mod.mu=as.numeric(coef(second.mod.est))[1]
    
    final.mod.coef=second.mod.coef
    nonzero.select=which(final.mod.coef!=0)
    final.mod.int=second.mod.mu
    final.mod.s2=sum((second.y-predict(second.mod.est,newx=second.modX))^2)/
                 (length(second.y)-sum(final.mod.coef[nonzero.select]!=0)-1)
    
    out=list(final.mod.coef=final.mod.coef, #Final Selection of Coefficients
            final.mod.int=final.mod.int,   #Final Estimated Intercept
            final.mod.s2=final.mod.s2,     #Final Estimated Noise Variance
            nonzero.select=nonzero.select) #Identifies the Nonzero Parameters
  }

  return(out)
}

#Illustration of Function for ADLASSO Estimation
adlasso6=adshrink678.func(x=maunaloa.co2.train,h=1,
         long.ar.select=F,maxP=14,maxQ=14,
         updateMA=F,BIC1=F,BIC2=F,eta=2,alpha=seq(0,1,0.1),
         Method="ADLASSO",K=5)
adlasso7=adshrink678.func(x=maunaloa.co2.train,h=1,
         long.ar.select=F,maxP=14,maxQ=14,
         updateMA=F,BIC1=F,BIC2=T,eta=2,alpha=seq(0,1,0.1),
         Method="ADLASSO",K=10)
adlasso8=adshrink678.func(x=maunaloa.co2.train,h=1,
         long.ar.select=F,maxP=14,maxQ=14,
         updateMA=F,BIC1=F,BIC2=T,eta=2,alpha=seq(0,1,0.1),
         Method="ADLASSO")

#Illustration of Function for ADENET Estimation
adenet6=adshrink678.func(x=maunaloa.co2.train,h=1,
        long.ar.select=F,maxP=14,maxQ=14,
        updateMA=F,BIC1=F,BIC2=F,eta=2,alpha=seq(0,1,0.1),
        ADENET.final="min",Method="ADENET",K=5)
adenet7=adshrink678.func(x=maunaloa.co2.train,h=1,
                        long.ar.select=F,maxP=14,maxQ=14,
        updateMA=F,BIC1=F,BIC2=T,eta=2,alpha=seq(0,1,0.1),
        ADENET.final="min",Method="ADENET",K=10)
adenet8=adshrink678.func(x=maunaloa.co2.train,h=1,
        long.ar.select=F,maxP=14,maxQ=14,
        updateMA=F,BIC1=F,BIC2=T,eta=2,alpha=seq(0,1,0.1),
        ADENET.final="min",Method="ADENET")
###############################################################################



###############################################################################
#Functions Used  to Perform K-fold Non-Dependent Cross Validation. 
#   All Three methods are based on dividing time series into blocks and 
#   performing CV with the blocks removing data that has mutual dependence
#   with test and training sets
#Arguments: x = time series
#           max.pq = maximum ar/ma order to consider 
#           (covers temporal dependence)
#           K = Number of Folds to Consider
#Source: Burman(2000), Racine(2000), and Bergmeir(2018)
###############################################################################


#Method 1: My Method 1 (Partition Data Into Blocks of Length max.pq)
NonDepCV1.func<-function(x,max.pq=NULL,K=NULL){
  N=length(x)
  if(is.null(max.pq)) max.pq=ceiling(10*log10(N))
  if(max.pq>ceiling(10*log10(N))){
    warning("Maximum ARMA Order Considered 
            Too Large Based on Time Series Length")
  }
  nblocks=floor(N/max.pq)
  if(is.null(K)) K=nblocks
  blocks=rep(nblocks,N)
  for(b in 1:(nblocks-1)){
    blocks[(b*max.pq-max.pq+1):(b*max.pq)]=b
  }
  if(K>nblocks){ 
    warning("Maximum Number of Folds Surpassed Based on Choice of max.pq
            and Time Series Length")
  }
  test.blocks=sort(sample(1:nblocks,K))
  Block.Matrix=matrix(0,N,K) 
  for(v in 1:K){
    in.test=which(blocks==test.blocks[v]) 
    Block.Matrix[in.test,v]=1 
    if(test.blocks[v]==1){
      out.ignore=which(blocks==2)
      Block.Matrix[out.ignore,v]=NA
    }else if(test.blocks[v]==nblocks){
      out.ignore=which(blocks==(nblocks-1))
      Block.Matrix[out.ignore,v]=NA
    }else{
      out.ignore=which(blocks %in% c(test.blocks[v]-1,test.blocks[v]+1))
      Block.Matrix[out.ignore,v]=NA
    }
  }
  return(Block.Matrix)
} 

#Method 2: Begmeir Method (Divide Data Into K Blocks and 
#          Remove Boundary Dependent Data)
NonDepCV2.func<-function(x,max.pq=NULL,K=NULL){
  N=length(x)
  obs.per.block=floor(N/K)
  if(is.null(max.pq)) max.pq=ceiling(10*log10(obs.per.block))
  if((max.pq/2)>obs.per.block){
    warning("Choice of max.pq and K Lead to No Observations")
  }
  blocks=rep(NA,N)
  for(b in 1:K){
    blocks[(b*obs.per.block-obs.per.block+
              ceiling(max.pq/2)+1):(b*obs.per.block-ceiling(max.pq)/2)]=b
  }
  
  Block.Matrix=matrix(NA,N,K) 
  for(v in 1:K){
    in.test=which(blocks==(1:K)[v]) 
    in.train=which(blocks %in% (1:K)[-v])
    Block.Matrix[in.test,v]=1 
    Block.Matrix[in.train,v]=0
  }
  
  return(Block.Matrix)
}    
#############################################################################



#############################################################################
#Function to Perform ADLASSO or ADENET Subse ARMA Estimation With
#Tuning Parameter Selection Based on Dependent K-fold Cross Validation 
#
#Arguments: x = Time Series to Be Modeled Using subset ARMA(maxP,maxQ)
#           h = Horizon Specific Model (Defaults to 1)
#           long.ar.select = Indicator Determining if Model Selection Is 
#                            Performed in the Initial Modeling of 
#                           Long AR Process(Defaults to F)
#           maxP = Maximum AR Order 
#           maxQ = Maximum MA Order
#           updateMA = Indicator Determining if Moving Average Terms Is 
#                      Updated After Initial Coefficients Selected 
#                      (Defaults to F)
#           alpha = Elastic Net Mixing Parameter 
#                   (0=Ridge,1=Lasso,Other=Elastic Net)
#                   (Defaults to Sequence 0,0.1,0.2,...,1)
#           eta = Exponent Applied to Weights (Defaults to 2)
#           Method = Choose Either "ADLASSO" or "ADENET"
#           K = Number of Folds for General CV (Defaults to NULL => LOOCV)
#           max.pq = Maximimum Temporal Dependence Considered in depOOS
#                    (Defaults to max(maxP,maxQ))
#           CV = Choose Either "KFOLD" or "LOBOCV"
#           ADENET.final = Appoach for Choosing Final Lambda for Each Alpha
#                          "min" = Choose Based on Minimum MSE
#                          "1se" = Choose Largest Lambda that Results in MSE
#                                  within 1 Standard Error of the Minimum
#
#Source: Wang and Leng(2007) and Efron et al.(2004) and Zou and Hastie(2005)
##############################################################################

#Creation of Function
adshrink91011.func<-function(x,h=1,long.ar.select=F,maxP,maxQ,updateMA=F,
                        eta=2,alpha=seq(0,1,0.1),
                        ADENET.final=c("min","1se"),max.pq = max(maxP,maxQ),
                        K=NULL,Method=c("ADLASSO","ADENET"),
                        CV=c("KFOLD","LOBOCV")){
  
  #Package Required
  require(glmnet) #Performs Ridge, Lasso, Elastic Net Estimation
  
  Method=match.arg(Method)
  CV=match.arg(CV)
  ADENET.final=match.arg(ADENET.final)
  
  Nt=length(x) #Length of Input Time Series
  
  #Fit Long AR Model to Estimate Innovations
  max.ar.order=ceiling(10*log10(Nt)) #Maximum Autoregressive Order
  init.mod.est=ar(x,aic=long.ar.select, #Allows Stepwise Selection
                  order.max=max.ar.order,demean=T)
  init.mod.error=residuals(init.mod.est)
  init.mod.order=length(which(is.na(init.mod.error)))
  
  #Create Model Matrix of AR and MA terms
  dataP=foreach(p=1:maxP,.combine=cbind)%do%{
    lag.func(x,k=(p+h-1))
  }
  dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
    lag.func(init.mod.error,k=(q+h-1))
  }
  first.modX=as.matrix(cbind(dataP,
             dataQ))[-(1:(init.mod.order+max(maxP,maxQ)+h-1)),]
  first.y=x[-(1:(init.mod.order+max(maxP,maxQ)+h-1))]
  
  #Number of Alphas For Elastic Net
  n.alpha=length(alpha)
  
  #Estimation Via ADLASSO
  if(Method=="ADLASSO"){
    if(CV=="KFOLD"){
      Fold.Matrix=NonDepCV2.func(x=first.y,max.pq=max.pq,K=K)
      nfolds=dim(Fold.Matrix)[2]
      
      lambda1.seq=cv.glmnet(y=first.y,x=first.modX,
                            standardize=T,alpha=1)$lambda
      
      SQDEV1=foreach(f=1:nfolds,.combine=rbind)%do%{
        in.train=which(Fold.Matrix[,f]==0)
        in.test=which(Fold.Matrix[,f]==1)
        first.mod.est=glmnet(y=first.y[in.train],x=first.modX[in.train,],
                      standardize=T,alpha=1,lambda=lambda1.seq)
        first.mod.res=(first.y[in.test]-predict(first.mod.est,
                       newx=first.modX[in.test,]))^2
        first.mod.res
      }
      
      CVM1=apply(SQDEV1,2,mean)
      lambda1.min=lambda1.seq[which.min(CVM1)]
      lambda1.1se=lambda1.seq[min(which(CVM1<(min(CVM1)+
                  sd(CVM1)/sqrt(length(CVM1)))))]
      
      first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,
                    alpha=1,lambda=lambda1.min)
      first.mod.coef=as.numeric(coef(first.mod.est))[-1]
      first.mod.mu=as.numeric(coef(first.mod.est))[1]
      weights=abs(first.mod.coef+1/length(first.y))^(-eta)
      
      if(updateMA){
        update.mod.predict=rep(NA,length(x))
        update.mod.error=rep(0,length(x))
        for(v in (h+max(maxP,maxQ)):Nt){
          update.mod.predict[v]=first.mod.mu+
            x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
            update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
          update.mod.error[v]=x[v]-update.mod.predict[v]
        }
        update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
          lag.func(update.mod.error,k=(q+h-1))
        }
        second.modX=as.matrix(cbind(dataP,
                    update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
        second.y=x[-(1:(max(maxP,maxQ)+h-1))]
      }else{
        second.modX=first.modX
        second.y=first.y
      }
      
      lambda2.seq=cv.glmnet(y=second.y,x=second.modX,standardize=T,
                  parallel=F, alpha=1,penalty.factor=weights)$lambda
      
      SQDEV2=foreach(f=1:nfolds,.combine=rbind)%do%{
        in.train=which(Fold.Matrix[,f]==0)
        in.test=which(Fold.Matrix[,f]==1)
        first.mod.est=glmnet(y=second.y[in.train],x=second.modX[in.train,],
                      standardize=T,alpha=1,lambda=lambda2.seq,
                      penalty.factor=weights)
        first.mod.res=(second.y[in.test]-predict(first.mod.est,
                      newx=second.modX[in.test,]))^2
        first.mod.res
      }
      
      CVM2=apply(SQDEV2,2,mean)
      lambda2.min=lambda2.seq[which.min(CVM2)]
      lambda2.1se=lambda2.seq[min(which(CVM2<(min(CVM2)+
                  sd(CVM2)/sqrt(length(CVM2)))))]
      
      final.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,
                    alpha=1,lambda=lambda2.1se,penalty.factor=weights)
      final.mod.coef=as.numeric(coef(final.mod.est))[-1]
      nonzero.select=which(final.mod.coef!=0)
      final.mod.int=as.numeric(coef(final.mod.est))[1]
      final.mod.s2=sum((second.y-predict(final.mod.est,
                   newx=second.modX,s=lambda2.1se,
                   method="lambda"))^2)/(length(second.y)-
                   sum(final.mod.coef[nonzero.select]!=0)-1)
      
      out=list(final.mod.coef=final.mod.coef, #Final Selection of Coefficients
               final.mod.int=final.mod.int,   #Final Estimated Intercept
               final.mod.s2=final.mod.s2,     #Final Estimated Noise Variance
               nonzero.select=nonzero.select) #Identifies Nonzero Parameters
    }else{
      Fold.Matrix=NonDepCV1.func(first.y,max.pq=max.pq,K=K)
      nfolds=dim(Fold.Matrix)[2]
      
      lambda1.seq=cv.glmnet(y=first.y,x=first.modX,
                            standardize=T,alpha=1)$lambda
      
      SQDEV1=foreach(f=1:nfolds,.combine=rbind)%do%{
        in.train=which(Fold.Matrix[,f]==0)
        in.test=which(Fold.Matrix[,f]==1)
        first.mod.est=glmnet(y=first.y[in.train],x=first.modX[in.train,],
                      standardize=T,alpha=1,lambda=lambda1.seq)
        first.mod.res=(first.y[in.test]-predict(first.mod.est,
                       newx=first.modX[in.test,]))^2
        first.mod.res
      }
      
      CVM1=apply(SQDEV1,2,mean)
      lambda1.min=lambda1.seq[which.min(CVM1)]
      lambda1.1se=lambda1.seq[min(which(CVM1<(min(CVM1)+
                  sd(CVM1)/sqrt(length(CVM1)))))]
      
      first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,
                           alpha=1,lambda=lambda1.min)
      first.mod.coef=as.numeric(coef(first.mod.est))[-1]
      first.mod.mu=as.numeric(coef(first.mod.est))[1]
      weights=abs(first.mod.coef+1/length(first.y))^(-eta)
      
      update.mod.predict=rep(NA,length(x))
      update.mod.error=rep(0,length(x))
      for(v in (h+max(maxP,maxQ)):Nt){
        update.mod.predict[v]=first.mod.mu+
          x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
          update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
        update.mod.error[v]=x[v]-update.mod.predict[v]
      }
      update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
        lag.func(update.mod.error,k=(q+h-1))
      }
      
      if(updateMA){
        update.mod.predict=rep(NA,length(x))
        update.mod.error=rep(0,length(x))
        for(v in (h+max(maxP,maxQ)):Nt){
          update.mod.predict[v]=first.mod.mu+
            x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
            update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
          update.mod.error[v]=x[v]-update.mod.predict[v]
        }
        update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
          lag.func(update.mod.error,k=(q+h-1))
        }
        second.modX=as.matrix(cbind(dataP,
                    update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
        second.y=x[-(1:(max(maxP,maxQ)+h-1))]
      }else{
        second.modX=first.modX
        second.y=first.y
      }
      
      lambda2.seq=cv.glmnet(y=second.y,x=second.modX,standardize=T,parallel=F,
                  alpha=1,penalty.factor=weights)$lambda
      SQDEV2=foreach(f=1:nfolds,.combine=rbind)%do%{
        in.train=which(Fold.Matrix[,f]==0)
        in.test=which(Fold.Matrix[,f]==1)
        first.mod.est=glmnet(y=second.y[in.train],x=second.modX[in.train,],
                      standardize=T,alpha=1,lambda=lambda2.seq,
                      penalty.factor=weights)
        first.mod.res=(second.y[in.test]-predict(first.mod.est,
                      newx=second.modX[in.test,]))^2
        first.mod.res
      }
      
      CVM2=apply(SQDEV2,2,mean)
      lambda2.min=lambda2.seq[which.min(CVM2)]
      lambda2.1se=lambda2.seq[min(which(CVM2<(min(CVM2)+
                  sd(CVM2)/sqrt(length(CVM2)))))]
      
      final.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,
                    alpha=1,lambda=lambda2.1se,penalty.factor=weights)
      final.mod.coef=as.numeric(coef(final.mod.est))[-1]
      nonzero.select=which(final.mod.coef!=0)
      final.mod.int=as.numeric(coef(final.mod.est))[1]
      final.mod.s2=sum((second.y-predict(final.mod.est,
                   newx=second.modX,s=lambda2.1se,
                   method="lambda"))^2)/(length(second.y)-
                   sum(final.mod.coef[nonzero.select]!=0)-1)
      
      out=list(final.mod.coef=final.mod.coef, #Final Selection of Coefficients
               final.mod.int=final.mod.int,   #Final Estimated Intercept
               final.mod.s2=final.mod.s2,     #Final Estimated Noise Variance
               nonzero.select=nonzero.select) #Identifies Nonzero Parameters
    }
  }
  
  #Estimation Via ADENET
  if(Method=="ADENET"){
    if(CV=="KFOLD"){
      Fold.Matrix=NonDepCV2.func(first.y,max.pq=max.pq,K=K)
      nfolds=dim(Fold.Matrix)[2]
      
      first.cv.out=NULL
      lambda1.seq=cv.glmnet(parallel=F,y=first.y,x=first.modX,standardize=T,
                  alpha=1)$lambda
      SQDEV1=foreach(f=1:nfolds,.combine=rbind)%do%{
        in.train=which(Fold.Matrix[,f]==0)
        in.test=which(Fold.Matrix[,f]==1)
        first.mod.est=glmnet(y=first.y[in.train],x=first.modX[in.train,],
                      standardize=T,alpha=1,lambda=lambda1.seq)
        first.mod.res=(first.y[in.test]-predict(first.mod.est,
                      newx=first.modX[in.test,]))^2
        first.mod.res
      }
      CVM1=apply(SQDEV1,2,mean)
      lambda1.min=lambda1.seq[which.min(CVM1)]
      lambda1.1se=lambda1.seq[min(which(CVM1<(min(CVM1)+
                  sd(CVM1)/sqrt(length(CVM1)))))]
      first.cv.out=rbind(first.cv.out,c(1,lambda1.min,min(CVM1)))
      
      
      first.mod.alpha=1
      first.mod.lambda=first.cv.out[which.min(first.cv.out[,3]),2]
      first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,
                    alpha=first.mod.alpha,lambda=first.mod.lambda)
      first.mod.coef=as.numeric(coef(first.mod.est))[-1]
      first.mod.mu=as.numeric(coef(first.mod.est))[1]
      
      weights=abs(first.mod.coef+1/length(first.y))^(-eta)
      
      if(updateMA){
        update.mod.predict=rep(NA,length(x))
        update.mod.error=rep(0,length(x))
        for(v in (h+max(maxP,maxQ)):Nt){
          update.mod.predict[v]=first.mod.mu+
            x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
            update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
          update.mod.error[v]=x[v]-update.mod.predict[v]
        }
        update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
          lag.func(update.mod.error,k=(q+h-1))
        }
        second.modX=as.matrix(cbind(dataP,
                    update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
        second.y=x[-(1:(max(maxP,maxQ)+h-1))]
      }else{
        second.modX=first.modX
        second.y=first.y
      }
      
      second.cv.out=NULL
      for(a in 1:n.alpha){
        lambda2.seq=cv.glmnet(parallel=F,y=second.y,x=second.modX,standardize=T,
                    alpha=alpha[a],penalty.factor=weights)$lambda
        SQDEV2=foreach(f=1:nfolds,.combine=rbind)%do%{
          in.train=which(Fold.Matrix[,f]==0)
          in.test=which(Fold.Matrix[,f]==1)
          second.mod.est=glmnet(y=second.y[in.train],x=second.modX[in.train,],
                  standardize=T,alpha=alpha[a],
                  penalty.factor=weights,lambda=lambda2.seq)
          second.mod.res=(second.y[in.test]-predict(second.mod.est,
                         newx=second.modX[in.test,]))^2
          second.mod.res
        }
        CVM2=apply(SQDEV2,2,mean)
        lambda2.min=lambda2.seq[which.min(CVM2)]
        lambda2.1se=lambda2.seq[min(which(CVM2<(min(CVM2)+
                    sd(CVM2)/sqrt(length(CVM2)))))]
        if(ADENET.final=="min") second.cv.out=rbind(second.cv.out,c(alpha[a],
                       lambda2.min,min(CVM2)))
        if(ADENET.final=="1se") second.cv.out=rbind(second.cv.out,c(alpha[a],
                       lambda2.1se,CVM2[min(which(CVM2<(min(CVM2)+
                       sd(CVM2)/sqrt(length(CVM2)))))]))
      }
      
      second.mod.alpha=alpha[which.min(second.cv.out[,3])]
      second.mod.lambda=second.cv.out[which.min(second.cv.out[,3]),2]
      second.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,
                     alpha=second.mod.alpha,lambda=second.mod.lambda,
                     penalty.factor=weights)
      second.mod.coef=as.numeric(coef(second.mod.est))[-1]
      second.mod.mu=as.numeric(coef(second.mod.est))[1]
      
      final.mod.coef=second.mod.coef
      nonzero.select=which(final.mod.coef!=0)
      final.mod.int=second.mod.mu
      final.mod.s2=sum((second.y-predict(second.mod.est,newx=second.modX))^2)/
                   (length(second.y)-sum(final.mod.coef[nonzero.select]!=0)-1)
      
      out=list(final.mod.coef=final.mod.coef, #Final Selection of Coefficients
               final.mod.int=final.mod.int,   #Final Estimated Intercept
               final.mod.s2=final.mod.s2,     #Final Estimated Noise Variance
               nonzero.select=nonzero.select) #Identifies Nonzero Parameters
    }else{
      Fold.Matrix=NonDepCV1.func(first.y,max.pq=max.pq,K=K)
      nfolds=dim(Fold.Matrix)[2]
      
      first.cv.out=NULL
      lambda1.seq=cv.glmnet(parallel=F,y=first.y,x=first.modX,standardize=T,
                  alpha=1)$lambda
      SQDEV1=foreach(f=1:nfolds,.combine=rbind)%do%{
        in.train=which(Fold.Matrix[,f]==0)
        in.test=which(Fold.Matrix[,f]==1)
        first.mod.est=glmnet(y=first.y[in.train],x=first.modX[in.train,],
                      standardize=T,alpha=1,lambda=lambda1.seq)
        first.mod.res=(first.y[in.test]-predict(first.mod.est,
                      newx=first.modX[in.test,]))^2
        first.mod.res
      }
      CVM1=apply(SQDEV1,2,mean)
      lambda1.min=lambda1.seq[which.min(CVM1)]
      lambda1.1se=lambda1.seq[min(which(CVM1<(min(CVM1)+
                  sd(CVM1)/sqrt(length(CVM1)))))]
      first.cv.out=rbind(first.cv.out,c(1,lambda1.min,min(CVM1)))
      
      
      first.mod.alpha=1
      first.mod.lambda=first.cv.out[which.min(first.cv.out[,3]),2]
      first.mod.est=glmnet(y=first.y,x=first.modX,standardize=T,
                    alpha=first.mod.alpha,lambda=first.mod.lambda)
      first.mod.coef=as.numeric(coef(first.mod.est))[-1]
      first.mod.mu=as.numeric(coef(first.mod.est))[1]
      
      weights=abs(first.mod.coef+1/length(first.y))^(-eta)
      
      if(updateMA){
        update.mod.predict=rep(NA,length(x))
        update.mod.error=rep(0,length(x))
        for(v in (h+max(maxP,maxQ)):Nt){
          update.mod.predict[v]=first.mod.mu+
            x[(v-h):(v-maxP-h+1)]%*%first.mod.coef[1:maxP]+
            update.mod.error[(v-h):(v-maxQ-h+1)]%*%first.mod.coef[-(1:maxP)]
          update.mod.error[v]=x[v]-update.mod.predict[v]
        }
        update.dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
          lag.func(update.mod.error,k=(q+h-1))
        }
        second.modX=as.matrix(cbind(dataP,
                    update.dataQ))[-(1:(max(maxP,maxQ)+h-1)),]
        second.y=x[-(1:(max(maxP,maxQ)+h-1))]
      }else{
        second.modX=first.modX
        second.y=first.y
      }
      
      second.cv.out=NULL
      for(a in 1:n.alpha){
        lambda2.seq=cv.glmnet(parallel=F,y=second.y,x=second.modX,standardize=T,
                    alpha=alpha[a],penalty.factor=weights)$lambda
        SQDEV2=foreach(f=1:nfolds,.combine=rbind)%do%{
          in.train=which(Fold.Matrix[,f]==0)
          in.test=which(Fold.Matrix[,f]==1)
          second.mod.est=glmnet(y=second.y[in.train],x=second.modX[in.train,],
                         standardize=T,alpha=alpha[a],penalty.factor=weights,
                         lambda=lambda2.seq)
          second.mod.res=(second.y[in.test]-predict(second.mod.est,
                          newx=second.modX[in.test,]))^2
          second.mod.res
        }
        CVM2=apply(SQDEV2,2,mean)
        lambda2.min=lambda2.seq[which.min(CVM2)]
        lambda2.1se=lambda2.seq[min(which(CVM2<(min(CVM2)+
                    sd(CVM2)/sqrt(length(CVM2)))))]
        if(ADENET.final=="min") second.cv.out=rbind(second.cv.out,
                                c(alpha[a],lambda2.min,min(CVM2)))
        if(ADENET.final=="1se") second.cv.out=rbind(second.cv.out,c(alpha[a],
                                lambda2.1se,CVM2[min(which(CVM2<(min(CVM2)+
                                sd(CVM2)/sqrt(length(CVM2)))))]))
      }
      
      second.mod.alpha=alpha[which.min(second.cv.out[,3])]
      second.mod.lambda=second.cv.out[which.min(second.cv.out[,3]),2]
      second.mod.est=glmnet(y=second.y,x=second.modX,standardize=T,
                     alpha=second.mod.alpha,lambda=second.mod.lambda,
                     penalty.factor=weights)
      second.mod.coef=as.numeric(coef(second.mod.est))[-1]
      second.mod.mu=as.numeric(coef(second.mod.est))[1]
      
      final.mod.coef=second.mod.coef
      nonzero.select=which(final.mod.coef!=0)
      final.mod.int=second.mod.mu
      final.mod.s2=sum((second.y-predict(second.mod.est,newx=second.modX))^2)/
        (length(second.y)-sum(final.mod.coef[nonzero.select]!=0)-1)
      
      out=list(final.mod.coef=final.mod.coef, #Final Selection of Coefficients
               final.mod.int=final.mod.int,   #Final Estimated Intercept
               final.mod.s2=final.mod.s2,     #Final Estimated Noise Variance
               nonzero.select=nonzero.select) #Identifies Nonzero Parameters
    }
  }
  
  return(out)
}

#Illustration of Function for ADLASSO Estimation
adlasso9=adshrink91011.func(x=maunaloa.co2.train,h=1,
         long.ar.select=F,maxP=14,maxQ=14,
         updateMA=F,eta=2,alpha=seq(0,1,0.1),max.pq=max(14,14),
         Method="ADLASSO",K=5,CV="KFOLD")
adlasso10=adshrink91011.func(x=maunaloa.co2.train,h=1,
          long.ar.select=F,maxP=14,maxQ=14,
          updateMA=F,eta=2,alpha=seq(0,1,0.1),max.pq=max(14,14),
          Method="ADLASSO",K=10,CV="KFOLD")
adlasso11=adshrink91011.func(x=maunaloa.co2.train,h=1,
          long.ar.select=F,maxP=14,maxQ=14,
          updateMA=F,eta=2,alpha=seq(0,1,0.1),max.pq=max(14,14),
          Method="ADLASSO",CV="LOBOCV")

#Illustration of Function for ADENET Estimation
adenet9=adshrink91011.func(x=maunaloa.co2.train,h=1,
        long.ar.select=F,maxP=14,maxQ=14,
        updateMA=F,eta=2,alpha=seq(0,1,0.1),max.pq=max(14,14),
        ADENET.final="min",Method="ADENET",K=5,CV="KFOLD")
adenet10=adshrink91011.func(x=maunaloa.co2.train,h=1,
         long.ar.select=F,maxP=14,maxQ=14,
         updateMA=F,eta=2,alpha=seq(0,1,0.1),max.pq=max(14,14),
         ADENET.final="min",Method="ADENET",K=10,CV="KFOLD")
adenet11=adshrink91011.func(x=maunaloa.co2.train,h=1,
         long.ar.select=F,maxP=14,maxQ=14,
         updateMA=F,eta=2,alpha=seq(0,1,0.1),max.pq=max(14,14),
         ADENET.final="min",Method="ADENET",CV="LOBOCV")
##############################################################################



########################################################################
#Function to Obtain Projected Posterior Distribution from Full Posterior
#Arguments: fullpost = Posterior Distribution from BHS estimated model
#           X = Model Matrix of Full Model
#           indproj = Vector Indicating which Columns are Included
#Source: Piironen and Vehtari (2015)
########################################################################

#Essential Function for Predictive Posterior Projection Method
pms.proj.func<-function(fullpost,X,indproj){
  
  #Transpose Matrix of Posterior Samples of ARMA Coefficients
  lres.coef=t(fullpost[,-dim(fullpost)[2]]) 
  
  #Vector of Posterior Samples of Variance Parameter
  lres.s2=(fullpost[,dim(fullpost)[2]])  
  S=length(lres.s2)   #Number of Posterior Samples Obtained
  N=dim(X)[1]   #Number of Points in Data Set
  
  P=dim(X)[2]
  COEF.PROJ=matrix(0,S,P)
  X.proj=X[,indproj]
  pred.proj=X%*%lres.coef
  
  coef.proj=solve(t(X.proj)%*%X.proj)%*%t(X.proj)%*%pred.proj
  var.proj=c(lres.s2)+colMeans((pred.proj-X.proj%*%coef.proj)^2)
  KL.PROJ=0.5*log(var.proj/lres.s2)
  COEF.PROJ[,indproj]=t(coef.proj)
  KL.MEAN=mean(KL.PROJ)
  COEF.MEAN=colMeans(COEF.PROJ)
  VAR.PROJ=var.proj
  VAR.MEAN=mean(var.proj)
  return(list(KL.MEAN=KL.MEAN,     #Average KL Divergence
              COEF.MEAN=COEF.MEAN, #Posterior Mean of Coefficients
              COEF.PROJ=COEF.PROJ, #Projected Posterior of Coefficients
              VAR.MEAN=VAR.MEAN,   #Posterior Mean of Variance
              VAR.PROJ=VAR.PROJ))  #Projected Posterior of Variance
} 
########################################################################



##############################################################################
#Function to Conduct ARMA Selection via Bayesian Projection Posterior 
#   Predictive Distribution Implementing Relative Efficiency
#   for Final Model Selection
#
#Arguments: x = Time Series to Be Modeled Using ARMA Process
#           h = Horizon Specific Model (Defaults to 1)
#           maxP = Maximum Autoregressive Order 
#           maxQ = Maximum Moving Average Order
#           updateMA = Indicator Determining if Moving Average Terms Should Be 
#                      Updated After Initial Coefficients Selected 
#                      (Defaults to F)
#           prior.choice = Choose Between Bayesian Horseshoe Prior ("hs") and 
#                          Bayesian Horseshoe+ Prior ("hs+")
#           KL.threshold = Single value or vector of chosen thresholds 
#                          for stopping rule based on Relative Efficiency 
#                          Comparing Submodel to Full Model
#                          based on Kullback Leibler Divergence
#                          (Defaults to c(0.9,0.95,0.99))
#Source: Piironen and Vehtari (2015)
############################################## ###############################

#Creation of Function
pms123.func<-function(x,h=1,maxP,maxQ,KL.threshold=c(0.90,0.95,0.99),
                       prior.choice=c("hs","hs+"),updateMA=F){
  require(MCMCpack)
  require(bayesreg)
  
  prior.choice=match.arg(prior.choice)
  
  N=length(x)
  max.ar.order=ceiling(10*log10(N))
  
  init.modX=foreach(init.ar=1:max.ar.order,.combine=cbind)%do%{
    lag.func(x,k=(init.ar+h-1))
  }

  init.data=data.frame(y=x,init.modX)
  init.data=init.data[-(1:(max.ar.order+h-1)),]
  
  init.mod.est=MCMCregress(y~.,data=init.data,mcmc=2000,thin=10,burnin=10000)
  muBeta0=mean(init.mod.est[,1])
  muBeta=colMeans(init.mod.est[,-c(1,dim(init.mod.est)[2])])
  init.mod.error=init.data$y-(as.numeric(muBeta0) + 
                 as.matrix(init.data[,-1])%*%as.vector(muBeta))
  
  dataP=foreach(p=1:maxP,.combine=cbind)%do%{
    lag.func(init.data$y,k=(p+h-1))
  }
  
  dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
    lag.func(init.mod.error,k=(q+h-1))
  }
  
  full.data=data.frame(y=init.data$y,dataP=dataP,dataQ=dataQ)
  full.data=full.data[-(1:(max(maxP,maxQ)+h-1)),]
  
  full.mod.est=bayesreg(y~.,data=full.data,prior=prior.choice,
               nsamples=2000,thin=10,burnin=10000)  
  full.mod.posterior=cbind(as.vector(full.mod.est$beta0),
                     t(as.matrix(full.mod.est$beta)),
                     as.vector(full.mod.est$sigma2))
  
  full.mod.int=c(full.mod.est$muBeta0)
  full.mod.coef=c(full.mod.est$muBeta)
  full.mod.s2=full.mod.est$muSigma2
  
  if(updateMA){
    full.mod.predict=rep(NA,length(x))
    full.mod.error=rep(0,length(x))
    for(k in (h+max(maxP,maxQ)):N){
      full.mod.predict[k]=full.mod.int+
        x[(k-h):(k-maxP-h+1)]%*%full.mod.coef[1:maxP]+
        full.mod.error[(k-h):(k-maxQ-h+1)]%*%full.mod.coef[-(1:maxP)]
      full.mod.error[k]=x[k]-full.mod.predict[k]
    }
  
    full.dataP=foreach(p=1:maxP,.combine=cbind)%do%{
      lag.func(x,k=(p+h-1))
    }
    
    update.dataQ=foreach(p=1:maxQ,.combine=cbind)%do%{
      lag.func(full.mod.error,k=(p+h-1))
    }
    
    full.mod.X=cbind(1,full.dataP,update.dataQ)
    full.mod.X=full.mod.X[-(1:(max(maxP,maxQ)+h-1)),]
    nMod=dim(full.mod.X)[2]
  }else{
    full.mod.X=as.matrix(cbind(1,full.data[,-1]))
    nMod=dim(full.mod.X)[2]
  }
  
  
  KL=rep(NA,nMod)
  chosen=1
  notchosen=setdiff(1:nMod,chosen)
  FIRST=pms.proj.func(fullpost=full.mod.posterior,X=full.mod.X,indproj=chosen)
  KL[1]=FIRST$KL.MEAN
  
  for(modnum in 2:nMod){
    nleft<-length(notchosen)
    val=foreach(j=1:nleft,.combine=c)%do%{
      ind<-sort(c(chosen,notchosen[j]))
      NEXT<-tryCatch({pms.proj.func(fullpost=full.mod.posterior,
                                    X=full.mod.X,indproj=ind)$KL.MEAN},
                     error=function(e){return(NA)})
      NEXT
    }
    minval<-which.min(val)
    chosen<-c(chosen,notchosen[minval])
    notchosen<-setdiff(1:nMod,chosen)
    KL[modnum]<-val[minval]
  }
  
  nKL.threshold=length(KL.threshold)
  KL.select=matrix(rep(chosen,nKL.threshold),ncol=nKL.threshold)
  final.mod.posterior=list()
  final.mod.int=list()
  final.mod.coef=list()
  final.mod.s2=list()
  
  for(k in 1:nKL.threshold){
    TEMP1=min(which((1-KL/KL[1])>KL.threshold[k]))
    TEMP2=pms.proj.func(fullpost=full.mod.posterior,
                        X=full.mod.X,indproj=chosen[1:TEMP1])
    KL.select[-(1:TEMP1),k]=NA
    final.mod.posterior[[k]]=cbind(TEMP2$COEF.PROJ,TEMP2$VAR.PROJ)
    final.mod.int[[k]]=TEMP2$COEF.MEAN[1]
    final.mod.coef[[k]]=TEMP2$COEF.MEAN[-1]
    final.mod.s2[[k]]=TEMP2$VAR.MEAN
  }
  
  return(list(CHOICE=chosen,  #Full Path of Variables in Order Of Selection
                              #  In Forward Search Algorithm
              KL=KL,  #Full Path of Kullback Leibler Divergences in Order
                      #  Of Selection in Forward Algorithm
              KL.threshold=KL.threshold, #Reoutputs the Thresholds Considered
              
              #Matrix where Each Column Identifies the Parameters Selected
              #Based on a Specific Threshold
              KL.select=KL.select, 
              
              #Estimated Intercept Before Selection                               
              full.mod.int=full.mod.int,
              #Estimated Coefficients Before Selection 
              full.mod.coef=full.mod.coef, 
              #Estimated Noise Variance Before Selection 
              full.mod.s2=full.mod.s2,       
              
              #Lists Where Each Element Corresponds to a
              #Specific Threshold in KL.threshold
              
              #Final Estimated Intercept for each Threshold
              final.mod.int=final.mod.int,
              #Final Selection of Coefficients for each Threshold
              final.mod.coef=final.mod.coef, 
              #Final Estimated Noise Variance for each Threshold
              final.mod.s2=final.mod.s2))    
}  

#Illustration of Function for BHS Estimation
bhs123.out=pms123.func(x=maunaloa.co2.train,h=1,maxP=14,maxQ=14,
                       KL.threshold=c(0.90,0.95,0.98),
                       prior.choice="hs",updateMA=F)
bhs1=list(final.mod.coef=bhs123.out$final.mod.coef[[1]], 
          final.mod.int=bhs123.out$final.mod.int[[1]],  
          final.mod.s2=bhs123.out$final.mod.s2[[1]],     
          nonzero.select=which(bhs123.out$final.mod.coef[[1]]!=0)) 
bhs2=list(final.mod.coef=bhs123.out$final.mod.coef[[2]], 
          final.mod.int=bhs123.out$final.mod.int[[2]],  
          final.mod.s2=bhs123.out$final.mod.s2[[2]],     
          nonzero.select=which(bhs123.out$final.mod.coef[[2]]!=0)) 
bhs3=list(final.mod.coef=bhs123.out$final.mod.coef[[3]], 
          final.mod.int=bhs123.out$final.mod.int[[3]],  
          final.mod.s2=bhs123.out$final.mod.s2[[3]],     
          nonzero.select=which(bhs123.out$final.mod.coef[[3]]!=0)) 

#Illustration of Function for BHS+ Estimation
bhsp123.out=pms123.func(x=maunaloa.co2.train,h=1,maxP=14,maxQ=14,
                        KL.threshold=c(0.90,0.95,0.98),
                        prior.choice="hs+",updateMA=F)
bhsp1=list(final.mod.coef=bhsp123.out$final.mod.coef[[1]], 
          final.mod.int=bhsp123.out$final.mod.int[[1]],  
          final.mod.s2=bhsp123.out$final.mod.s2[[1]],     
          nonzero.select=which(bhsp123.out$final.mod.coef[[1]]!=0)) 
bhsp2=list(final.mod.coef=bhsp123.out$final.mod.coef[[2]], 
          final.mod.int=bhsp123.out$final.mod.int[[2]],  
          final.mod.s2=bhsp123.out$final.mod.s2[[2]],     
          nonzero.select=which(bhsp123.out$final.mod.coef[[2]]!=0)) 
bhsp3=list(final.mod.coef=bhsp123.out$final.mod.coef[[3]], 
          final.mod.int=bhsp123.out$final.mod.int[[3]],  
          final.mod.s2=bhsp123.out$final.mod.s2[[3]],     
          nonzero.select=which(bhsp123.out$final.mod.coef[[3]]!=0)) 
#############################################################################



#############################################################################
#Function to Conduct ARMA Selection via 
#   Bayesian Projection Posterior Predictive Distribution
#   Implementing Out-of-Sample Based Final Model Selection
#
#Arguments: x = Time Series to Be Modeled Using ARMA Process
#           h = Horizon Specific Model (Defaults to 1)
#           maxP = Maximum Autoregressive Order 
#           maxQ = Maximum Moving Average Order
#           updateMA = Indicator Determining if Moving Average Terms Is 
#                      Updated After Initial Coefficients Selected 
#                      (Defaults to F)
#           KL.threshold = Single value or vector of chosen thresholds 
#                          for stopping rule based on Relative Efficiency 
#                          Comparing Submodel to Full Model
#                          based on Kullback Leibler Divergence
#                          (Defaults to c(0.9,0.95,0.99))
#           KL.stop = Stopping Rule 
#     (Select from Considered Models Where Relative Efficiency < KL.stop)
#Source: Piironen and Vehtari (2015)
#############################################################################

#Creation of Function
pms4.func<-function(x,h=1,maxP,maxQ,KL.stop=0.98,test.per=0.2,
                    prior.choice=c("hs","hs+"),updateMA=F){
  require(MCMCpack)
  require(bayesreg)
  
  prior.choice=match.arg(prior.choice)
  
  cv.vector=OOS.IndepCV.func(x,test.per=test.per)
  
  x.train=x[cv.vector==0]
  x.test=x[cv.vector==1]
  
  Nt=length(x.train)
  max.ar.order=ceiling(10*log10(Nt))
  
  init.modX=foreach(init.ar=1:max.ar.order,.combine=cbind)%do%{
    lag.func(x.train,k=(init.ar+h-1))
  }
  
  init.data=data.frame(y=x.train,init.modX)
  init.data=init.data[-(1:(max.ar.order+h-1)),]
  
  init.mod.est=MCMCregress(y~.,data=init.data,mcmc=2000,thin=10,burnin=10000)
  muBeta0=mean(init.mod.est[,1])
  muBeta=colMeans(init.mod.est[,-c(1,dim(init.mod.est)[2])])
  init.mod.error=init.data$y-(as.numeric(muBeta0) + 
                 as.matrix(init.data[,-1])%*%as.vector(muBeta))
  
  dataP=foreach(p=1:maxP,.combine=cbind)%do%{
    lag.func(init.data$y,k=(p+h-1))
  }
  
  dataQ=foreach(q=1:maxQ,.combine=cbind)%do%{
    lag.func(init.mod.error,k=(q+h-1))
  }
  
  full.data=data.frame(y=init.data$y,dataP=dataP,dataQ=dataQ)
  full.data=full.data[-(1:(max(maxP,maxQ)+h-1)),]
  
  #############################################
  #xc.mean=as.numeric(colMeans(full.data))
  #xc.sd=as.numeric(apply(full.data,2,sd))
  #full.data=as.data.frame(scale(full.data))
  ##############################################
  
  full.mod.est=bayesreg(y~.,data=full.data,prior=prior.choice,nsamples=2000,
                        thin=10,burnin=10000)  
  full.mod.posterior=cbind(as.vector(full.mod.est$beta0),
                           t(as.matrix(full.mod.est$beta)),
                           as.vector(full.mod.est$sigma2))
  
  full.mod.int=c(full.mod.est$muBeta0)
  full.mod.coef=c(full.mod.est$muBeta)
  full.mod.s2=full.mod.est$muSigma2
  
  if(updateMA){
    full.mod.predict=rep(NA,length(x.train))
    full.mod.error=rep(0,length(x.train))
    for(k in (h+max(maxP,maxQ)):Nt){
      full.mod.predict[k]=full.mod.int+
         x.train[(k-h):(k-maxP-h+1)]%*%full.mod.coef[1:maxP]+
         full.mod.error[(k-h):(k-maxQ-h+1)]%*%full.mod.coef[-(1:maxP)]
      full.mod.error[k]=x.train[k]-full.mod.predict[k]
    }
    
    full.dataP=foreach(p=1:maxP,.combine=cbind)%do%{
      lag.func(x.train,k=(p+h-1))
    }

    update.dataQ=foreach(p=1:maxQ,.combine=cbind)%do%{
      lag.func(full.mod.error,k=(p+h-1))
    }
    
    full.mod.X=cbind(1,full.dataP,update.dataQ)
    full.mod.X=full.mod.X[-(1:(max(maxP,maxQ)+h-1)),]
    nMod=dim(full.mod.X)[2]
  }else{
    full.mod.X=as.matrix(cbind(1,full.data[,-1]))
    nMod=dim(full.mod.X)[2]
  }
  
  KL=rep(NA,nMod)
  chosen=1
  notchosen=setdiff(1:nMod,chosen)
  FIRST=pms.proj.func(fullpost=full.mod.posterior,X=full.mod.X,indproj=chosen)
  KL[1]=FIRST$KL.MEAN
  modnum=2
  check=0
  
  while(check<KL.stop | modnum==nMod){
    nleft<-length(notchosen)

    val=foreach(j=1:nleft,.combine=c)%do%{
      ind<-sort(c(chosen,notchosen[j]))
      NEXT<-tryCatch({pms.proj.func(fullpost=full.mod.posterior,
                     X=full.mod.X,indproj=ind)$KL.MEAN},
                     error=function(e){return(NA)})
      NEXT
    }
    minval<-which.min(val)
    chosen<-c(chosen,notchosen[minval])
    notchosen<-setdiff(1:nMod,chosen)
    KL[modnum]<-val[minval]
    check=1-KL[modnum]/KL[1]
    modnum=modnum+1
  }
  
  KL=KL[!is.na(KL)]
  chosen=chosen[!is.na(chosen)]
  
  nMod2=length(KL)
  MSE=rep(NA,nMod2)
  
  for(modnum in 1:nMod2){
    out=pms.proj.func(fullpost=full.mod.posterior,X=full.mod.X,
                      indproj=chosen[1:modnum])
    coef=out$COEF.MEAN
    int=coef[1]
    coef.ar=coef[2:(1+maxP)]
    coef.ma=coef[-(1:(maxP+1))]
    
    predictx.test=rep(NA,length(x.test))
    errorx.test=rep(0,length(x.test))
    for(k in (h+max(maxP,maxQ)):length(x.test)){
      predictx.test[k]=int+x.test[(k-h):(k-maxP-h+1)]%*%coef.ar+
        errorx.test[(k-h):(k-maxQ-h+1)]%*%coef.ma
      errorx.test[k]=x.test[k]-predictx.test[k]
    }
    
    MSE[modnum]=mean((x.test-predictx.test)^2,na.rm=T)
  }
  best.mod=chosen[1:which.min(MSE)]
  out.mod=pms.proj.func(fullpost=full.mod.posterior,
                        X=full.mod.X,indproj=best.mod)
  
  final.mod.int=out.mod$COEF.MEAN[1]
  final.mod.coef=out.mod$COEF.MEAN[-1]
  final.mod.s2=out.mod$VAR.MEAN
  
  nonzero.select=which(final.mod.coef!=0)
  
  return(list(CHOICE=chosen,  #Full Path of Variables in Order Of Selection
              #  In Forward Search Algorithm
              KL=KL,  #Full Path of Kullback Leibler Divergences in Order
              #  Of Selection in Forward Algorithm
              KL.threshold=KL.threshold, #Reoutputs the Thresholds Considered
              
              #Matrix where Each Column Identifies the Parameters Selected
              #Based on a Specific Threshold
              KL.select=KL.select, 
              
              #Estimated Intercept Before Selection                               
              full.mod.int=full.mod.int,
              #Estimated Coefficients Before Selection 
              full.mod.coef=full.mod.coef, 
              #Estimated Noise Variance Before Selection 
              full.mod.s2=full.mod.s2,       
              
              #Lists Where Each Element Corresponds to a
              #Specific Threshold in KL.threshold
              
              #Final Estimated Intercept for each Threshold
              final.mod.int=final.mod.int,
              #Final Selection of Coefficients for each Threshold
              final.mod.coef=final.mod.coef, 
              #Final Estimated Noise Variance for each Threshold
              final.mod.s2=final.mod.s2))    
}   

#Illustration of Function for BHS Estimation
bhs4=pms4.func(x=maunaloa.co2.train,h=1,maxP=14,maxQ=14,KL.stop=0.98,
                          test.per=0.2,prior.choice="hs",updateMA=F)

#Illustration of Function for BHS+ Estimation
bhsp4=pms4.func(x=maunaloa.co2.train,h=1,maxP=14,maxQ=14,KL.stop=0.98,
                           test.per=0.2,prior.choice="hs+",updateMA=F)
#############################################################################



#############################################################################
#Obtain Forecasts for All Models
#
#List of All Models

MODELS=list(adlasso1,adlasso2,adlasso3,adlasso4,adlasso5,
            adlasso6,adlasso7,adlasso8,
            adlasso9,adlasso10,adlasso11,
            adenet1,adenet2,adenet3,adenet4,adenet5,
            adenet6,adenet7,adenet8,adenet9,
            adenet10,adenet11,
            bhs1,bhs2,bhs3,bhs4,
            bhsp1,bhsp2,bhsp3,bhsp4)

#
#############################################################################

#Total Number of Models Estimated
nMODELS=length(MODELS) 

#Each Row of the Following Matrices Corresponds to a Different Final Model

#Matrix of Binary Variables Indicated Selection
SELECT=matrix(0,nMODELS,ncol=28) 
#Matrix of Coefficients Estimated from All Models
COEF=matrix(0,nMODELS,ncol=28)

#1-step Ahead Forecasting Results
RMSFE=rep(NA,nMODELS) #Matrix of Root mean squared Forecast Error
MAPFE=rep(NA,nMODELS) #Matrix of Mean Absolute Percentage Forecast Error
MFB=rep(NA,nMODELS)   #Matrix of Mean Forecast Bias
MDFB=rep(NA,nMODELS)  #Matrix of Mean Directional Forecast Bias

#Matrices of Lower, Upper, and Mean Forecasts 
#   (Monte Carlo 1-step Ahead Forecast Distribution)
#   Columns Correspond for each of the models estimated

#5% Quantile of Monte Carlo Distribution
FC.LOWER=matrix(NA,sum(VALIDATION.PERIOD),nMODELS)
#Mean of Monte Carlo Distribution
FC.MEAN=matrix(NA,sum(VALIDATION.PERIOD),nMODELS)
#95% Quantile of Monte Carlo Distribution
FC.UPPER=matrix(NA,sum(VALIDATION.PERIOD),nMODELS) 

for(k in 1:nMODELS){
  temp.mod=MODELS[[k]]
  SELECT[k,temp.mod$nonzero.select]=1
  COEF[k,]=temp.mod$final.mod.coef
  
  fc.lower=rep(NA,length(maunaloa.co2.final))
  fc.mean=rep(NA,length(maunaloa.co2.final))
  fc.upper=rep(NA,length(maunaloa.co2.final))
  
  error=rep(0,length(maunaloa.co2.final))
  
  for(j in 30:length(maunaloa.co2.final)){
    fc=as.numeric(maunaloa.co2.final[(j-1):(j-14)]%*%
                    temp.mod$final.mod.coef[1:14]+
                    error[(j-1):(j-14)]%*%temp.mod$final.mod.coef[-(1:14)]) +
      rnorm(100000,mean=temp.mod$final.mod.int,sd=sqrt(temp.mod$final.mod.s2))
    fc.lower[j]=quantile(fc,0.05)
    fc.mean[j]=mean(fc)
    fc.upper[j]=quantile(fc,0.95)
    error[j]=maunaloa.co2.final[j]-fc.mean[j]
  }
  
  RMSFE[k]=sqrt(mean((maunaloa.co2.final[VALIDATION.PERIOD]-
           fc.mean[VALIDATION.PERIOD])^2))
  MAPFE[k]=100*mean(abs((maunaloa.co2.final[VALIDATION.PERIOD]-
           fc.mean[VALIDATION.PERIOD])/
           maunaloa.co2.final[VALIDATION.PERIOD]))
  MFB[k]=mean(maunaloa.co2.final[VALIDATION.PERIOD]-
         fc.mean[VALIDATION.PERIOD])
  MDFB[k]=(sum((maunaloa.co2.final[VALIDATION.PERIOD]-
           fc.mean[VALIDATION.PERIOD])>0)-
           sum((maunaloa.co2.final[VALIDATION.PERIOD]-
           fc.mean[VALIDATION.PERIOD])<0))/
           sum(VALIDATION.PERIOD)
  
  FC.LOWER[,k]=fc.lower[VALIDATION.PERIOD]
  FC.MEAN[,k]=fc.mean[VALIDATION.PERIOD]
  FC.UPPER[,k]=fc.upper[VALIDATION.PERIOD]
  
}
##############################################################################



##############################################################################
#Check Stationarity and Invertibility of Estimates From All Different Methods
#Based on the COEF matrix above
##############################################################################

#Check Stationarity and Invertibility of Estimates
stationarity=rep(NA,nMODELS)
invertibility=rep(NA,nMODELS)
for(k in 1:nMODELS){
  ar.coef=COEF[k,1:14]
  st.poly=c(1,-ar.coef)
  st.root=polyroot(st.poly)
  st.check=sum(abs(st.root)>1)==length(st.root)
  stationarity[k]=st.check
  
  ma.coef=COEF[k,-(1:14)]
  ma.poly=c(1,ma.coef)
  ma.root=polyroot(ma.poly)
  ma.check=sum(abs(ma.root)>1)==length(ma.root)
  invertibility[k]=ma.check
}

print(stationarity) #All True
print(invertibility) #Models (adlasso4,adlasso5,adenet4,adenet5) Fail
######################## #####################################################





















































