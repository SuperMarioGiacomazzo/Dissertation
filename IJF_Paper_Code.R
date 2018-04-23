###################################################################################
#Paper:"BAYESIAN SHRINKAGE ESTIMATES OF LOGISTIC SMOOTH TRANSITION AUTOREGRESSIONS"
#Authors:Mario Giacomazzo (Arizona State University) 
#        Yiannis Kamarianakis (Arizona State University)
#Year:2017
#
#Comments: This code requires a working installation of JAGS by Martyn Plummer
###################################################################################

####################
#Required R Packages
####################
library(runjags) #Needed for Calling JAGS through R
library(doParallel) #Needed for Parallelization of MCMC chains
library(datasets) #Contains Sunspot DATA
library(tsDyn) #Used for Frequentist Estimation of LSTAR Models

##########################################
#Simulate 100 Replicates of LSTAR(3) Model
##########################################

#Specify Autoregressive (AR) Coefficients in Low and High Regimes
##Regime Specific AR models must reflect Stationarity Within Regimes
a0<-0
a1<-0
a2<-0
a3<--0.6
b0<-0.02
b1<-0
b2<-0
b3<-0.75

#Specify Transition Function Parameters
slope<-120 
thresh<-0.02
delay<-2

#Specify Standard Deviation of Error Term
sigma<-0.02

#Simulation Information
S=3 #Number of Replications
N=1000 #Length of Each Simulated Time Series Desired
burn=2000 #Burn-in Size for Each Replication

#Function Used to Generate 1 Replication of LSTAR Model
generate.func<-function(x){
  set.seed(x) #Needed to Obtain Different Replications That Are Reproducible
  y=rnorm((N+burn),0,0.02) #Initialize Time Series
  e=rnorm((N+burn),0,sigma) #Create iid Errors
  for(i in 4:(N+burn)){
    wt=1/(1+exp(-slope*(y[i-delay]-thresh))) 
    y[i]=(a0+a1*y[i-1]+a2*y[i-2]+a3*y[i-3])*(1-wt)+
      (b0+b1*y[i-1]+b2*y[i-2]+b3*y[i-3])*wt+e[i]
  }
  return(y[-(1:burn)]) #Output Time Series After Beginning Burn-in Period
}

all.data<-lapply(1:S,generate.func) #Create List of Many Replications

###################################################
#Function Required for Obtaining Lagged Time Series
###################################################
lag.func<-function(x,k=1){
  t=length(x)
  y=c(rep(NA,t))
  for(i in (k+1):t){
    y[i]=x[i-k]
  }
  return(y)
}

##################################################
#Necessary Elements for MCMC Sampling for JAGS
#Used to Model LSTAR(4) Model 
#Using Bayesian Horseshoe Priors for AR Parameters
#With Dirichlet Prior for Threshold Variable
##################################################

#Function Used to Obtain Replication Specific Data
datafunction<-function(i){
  y=all.data[[i]]
  X=matrix(NA,nrow=length(y),ncol=4)
  for(j in 1:4){
    X[,j]=lag.func(y,k=j)
  }
  X=cbind(1,X)
  return(list(y=all.data[[i]], #The i-th Replicated Time Series
              #Model Matrix For Each Regime
              #(Made up of Lags 1 to 4 of Endogenous Series)
              X=X, 
              N=length(all.data[[i]]), #Length of Time Series
              #Minimum Hyperparameter for Uniform Prior on Threshold
              min.thresh=quantile(all.data[[i]],0.15),
              #Maximum Hyperparameter for Uniform Prior on Threshold
              max.thresh=quantile(all.data[[i]],0.85), 
              #Standard Deviation of Endogenous Time Series 
              #Used to Scale Slope for Threshold Variable
              sdy=sd(y),
              #Matrix Containing All Delays Considered for Threshold Variable
              X2=X[,-1],
              #Hyperparameter for Dirichlet Distribution
              #(length must equal number of columns in X2;
              #   elements must sum to 1)
              prop.prior=c(.25,.25,.25,.25))) 
                                              
}

#JAGS Model Represented as a String 
#(Horseshoe Priors are Used for Shrinkage and
#Dirichlet Used for Threshold Variable)
#Notice: We do not Monitor Tuning Parameters and
#We Monitor the Raw Unscaled Slope
MOD<-"model{
    #Likelihood Function (Starts at p+1 Which in Our Case is 5)
    for(i in 5:N){
y[i]~dnorm(mu[i],tau) 
w[i]<-1.0/(1.0+exp(-(preslope/sdy)*(inprod(prop[],X2[i,])-thresh)))
mu[i]<-(inprod(alpha[],X[i,]))*(1.0-w[i])+(inprod(beta[],X[i,]))*(w[i])
}

tau~dgamma(.001,.001) #Prior for Error Precision
preslope~dlnorm(3,1) #Prior for Scaled Transition Slope Parameter 
                     #(Can Use Gamma, Truncated Normal, etc.)
thresh~dunif(min.thresh,max.thresh) #Prior for Threshold Variable
#Prior for Weights of Linear Combination of Possible Threshold Variables
prop~ddirch(prop.prior) 

global.squared<-global^2 #Global Shrinkage Parameter

for(k in 1:5){
#Local Shrinkage Parameters for Low Regime
local1.squared[k]<-(local1[k])^2 
#Local Shrinkage Parameters for High Regime
local2.squared[k]<-(local2[k])^2 
#Priors for Shrinkage Parameters for Low Regime
alpha[k]~dnorm(0,tau/(global.squared*local1.squared[k])) 
#Priors for Shrinkage Parameters for High Regime
beta[k]~dnorm(0,tau/(global.squared*local2.squared[k])) 
}

#Bayesian Global-Local Prior Hierarchy Using Half-Cauchy Distributions
#t-Distribution with 1 df -> Cauchy
#T(0,) Truncates Distribution from 0 to Infinity
for(k in 1:5){
local1[k]~dt(0,1,1)T(0,)
local2[k]~dt(0,1,1)T(0,) 
}
global~dt(0,1,1)T(0,)

#Find Raw Unscaled Transition Slope Parameter
slope<-preslope/sdy

#modules# runjags
#monitor# tau,slope,thresh,alpha,beta,prop
}
"

#Function Used to Obtain Initial Values Needed for All Parameters
#Different Initial Values for Different Chains for AR Parameters
initsfunction<-function(chain){
  tau<-c(0.01,20,100)[chain]
  preslope<-c(20,50,100)[chain]
  thresh<-0.02
  prop<-c(0.25,.25,0.25,0.25)
  set.seed(chain)
  alpha<-rnorm(5,0,1)
  global<-0.5
  local1<-rep(0.5,5)
  local2<-rep(0.5,5)
  set.seed(chain+1)
  beta<-rnorm(5,0,1)
  .RNG.seed<-c(1,2,3)[chain]
  .RNG.name<-c("base::Super-Duper","base::Wichmann-Hill","base::Super-Duper")
  return(list(tau=tau,preslope=preslope,thresh=thresh,prop=prop,
              alpha=alpha,beta=beta,local1=local1,local2=local2,global=global,
              .RNG.seed=.RNG.seed,.RNG.name=.RNG.name))
}

##################################################
#MCMC Posterior Sampling for Each Replication
##################################################

#Parallelization is Used Across the Many Replications Using Foreach Package
cl2<-makeCluster(1) #Number of Clusters if Access to >3 Cores
registerDoParallel(cl2)

#Foreach Package Outputs as a List Where Each Element is a Different Replication
hs.out=foreach(v=1:S,.packages=c("runjags","parallel")) %dopar%{
  #Parallelization Used Also for Different MCMC Chains
  cl<-makeCluster(3)
  #Initialize JAGS Model for Specific Replication Using 3 Chains
  model2<-run.jags(MOD,data=datafunction(v),n.chains=3,inits=initsfunction,
            mutate=list(prec2sd, 'tau'),adapt=5000,burnin=10000,sample=1000,
            thin=10,method="rjparallel",method.options=list(cl=cl))
  #Obtain Initial Maximum PSRF Convergence Statistic 
  #for Chain Convergence Across All Parameters
  max.psrf=max(summary(model2)[,"psrf"],na.rm=T)
  #Obtain Initial Minimum Effective Sample Size Across All Saved Parameters
  min.ess=min(summary(model2)[,"SSeff"],na.rm=T) 
  i=1 #Identify This as Initialized Model
  
  #If Convergence is not Met, Then Update Model With 1000*i samples
  #Repeat Until Convergence is Met or 20 updates have occurred 
  #(Could Take a Long Time)
  #The Maximum Number of Updates is currently 20 but may be reduced
  while((max.psrf>1.05|min.ess<150)&i<20){
    model2<-extend.jags(model2,adapt=1000,burnin=0,
                        sample=1000*i,silent.jags=T)
    max.psrf=max(summary(model2)[,"psrf"],na.rm=T)
    min.ess=min(summary(model2)[,"SSeff"],na.rm=T)
    i=i+1
    print(max.psrf)
    
  }  
  stopCluster(cl)
  
  #For each replication we output a list containing the final model, 
  #convergence results and total computation time
  out=list(model2=model2,max.psrf=max.psrf,min.ess=min.ess,
           time=round(as.numeric(model2$timetaken)/60,1))
}
stopCluster(cl2)

##################################################
#Analyzing Output From Simulation Study
##################################################

#True Parameters for Simulated Nonlinear LSTAR(3) Under Assumption that p=4
true<-c(0,0,0,-0.6,0,0.02,0,0,0.75,0,0.02,120,0.02)
true.delay2<-c(0,1,0,0)

#Function to Output final PSRF Statistic Which Determines 
#if Convergence was Met for Each Replication
conv.func<-function(x){
  return(x$max.psrf)
}

#Check Convergence

#Obtain Max PSRF for Each Replication
hsdlp2.conv=unlist(lapply(hs.out,conv.func)) 
#Check Which Replications Converged
id.hsdlp2.conv=which(hsdlp2.conv<1.05)
#Calculate Convergence Percentage
hsdlp2.per=length(id.hsdlp2.conv)/100 


#Check # of Samples Required For Replications Where Convergence Was Met
hsdlp2.samples=rep(NA,100)
for (k in 1:100){ 
  #Obtain Number of Samples Required For Convergence For Each Replication
  hsdlp2.samples[k]=hs.out[[k]]$model2$sample 
}
#Replace Number of Samples with NA For Replications that Didn't Converge
hsdlp2.samples[-id.hsdlp2.conv]=NA 


#Get Tables of Posterior Estimates of Nonlinear Parameters
HSDLP2.EST.PARAMS=matrix(NA,100,13)
HSDLP2.EST.THVAR=matrix(NA, 100,4 )
for (k in 1:100){
  HSDLP2.EST.PARAMS[k,]=summary(hs.out[[k]]$model2)[c(4:13,18,2:3),"Mean"]
  HSDLP2.EST.THVAR[k,]=summary(hs.out[[k]]$model2)[14:17,"Mean"]
}

#Plotting the Threshold Variable for The Second Threshold Variable
png(file="hsthvar2.png",height=600,width=600)
z1=c(0,1,0,-1)
z2=c(1,0,-1,0)
par(mar=c(1.1,1.1,1.1,1.1))
plot(z1,z2,plot="n",pch=".",xlim=c(-1.2,1.2),ylim=c(-1.2,1.2),
     xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
points(x=1,y=0,pch=16,col="black",add=T)
polygon(z1,z2,border="black")
polygon(z1/2,z2/2,border="black")
polygon(z1/4,z2/4,border="black")
polygon(3*z1/4,3*z2/4,border="black")
text(0,1.1,expression(y[t-1]),cex=2,col="black")
text(1.18,0,expression(y[t-2]),cex=2,col="black")
text(0,-1.1,expression(y[t-3]),cex=2,col="black")
text(-1.16,0,expression(y[t-4]),cex=2,col="black")
text(-0.18,0.18,0.25,col="black")
text(-0.31,0.31,0.5,col="black")
text(-0.44,0.44,0.75,col="black")
text(-0.58,0.58,"1.00",col="black")

for(k in id.hsdlp2.conv){
  x=HSDLP2.EST.THVAR[k,]
  x1=x
  y1=x
  x1[1]=0
  y1[2]=0
  x1[3]=0            
  y1[3]=-x[3]
  x1[4]=-x[4]
  y1[4]=0
  polygon(x1,y1,col=rgb(176/255,48/255,96/255,0.1),border=NA)
}
dev.off()

#Obtain RMSE for Each Parameter Obtained 
#Using Posterior Estimates Across All Simulations
rmse.func<-function(x.est){
  x.sqdiff=(x.est-true)^2
  return(x.sqdiff)
}
HSDLP2.RMSE=sqrt(rowMeans(apply(HSDLP2.EST[id.hsdlp2.conv,],1,rmse.func)))

################################################
#Practical Application to Annual Sunspot Numbers
################################################
#Import Training and Testing Datasets of Sunspots
ss.year=window(sunspot.year,start=1700, end=1748)
ss.month1=window(aggregate(sunspot.month, FUN=mean),start=1749,end=1979)
#Create Training Dataset Using Years 1700-1979
Train=c(as.vector(ss.year),as.vector(ss.month1)) 
#Create Testing Dataset Using Years 1980-2006
Test=window(aggregate(sunspot.month, FUN=mean),start=1980,end=2006) 

#Apply Classical Square Root Transformation To Train and Test Data
Train.transform=2*(sqrt(Train+1)-1)
Test.transform=2*(sqrt(Test+1)-1)

##################################################
#Necessary Elements for MCMC Sampling With JAGS
#Specific To Sunspot Data Using Both
#Bayesian Lasso and Bayesian Horseshoe Priors
#Along with Dirichlet Prior for Threshold Variable
##################################################
#General Data Function Data Function
datafunction<-function(){
  y=Train.transform
  X=matrix(NA,nrow=length(y),ncol=10)
  for(j in 1:10){
    X[,j]=lag.func(y,k=j)
  }
  X=cbind(1,X)
  return(list(y=y,X=X,N=length(y),
              min.thresh=quantile(y,0.15),
              max.thresh=quantile(y,0.85),
              sdy=sd(y)*sqrt(length(y)-1)/sqrt(length(y))))
}

#Bayesian Lasso Model With Regime Specific Shrinkage Parameters 
MOD1<-"
model{
for(i in 11:N){
y[i]~dnorm(mu[i],tau)
w[i]<-1.0/(1.0+exp(-(preslope/sdy)*(y[i-2]-thresh)))
mu[i]<-(inprod(alpha[],X[i,]))*(1.0-w[i])+(inprod(beta[],X[i,]))*(w[i])
}

tau~dgamma(.001,.001)
preslope~dlnorm(3,1)
thresh~dunif(min.thresh,max.thresh)

for(k in 1:11){
alpha[k]~dnorm(0,tau*alphatau[k])
beta[k]~dnorm(0,tau*betatau[k])
}

lt1<-lambda1.squared/2
lt2<-lambda2.squared/2

for(k in 1:11){
alphatau[k]~dexp(lt1)
betatau[k]~dexp(lt2)
}
#Gamma Prior for Low Regime Shrinkage Parameter
lambda1.squared~dgamma(1,1.78) 
#Gamma Prior for High Regime Shrinkage Parameter
lambda2.squared~dgamma(1,1.78) 

#modules# runjags
#monitor# tau,preslope,thresh,alpha,beta
}
"
initsfunction1<-function(chain){
  tau<-c(0.01,20,100)[chain]
  preslope<-c(5,10,15)[chain]
  thresh<-10.8
  set.seed(chain)
  alpha<-rnorm(11,0,1)
  set.seed(chain+1)
  beta<-rnorm(11,0,1)
  alphatau<-rgamma(11,.01,.01)
  betatau<-rgamma(11,.01,.01)
  lambda1.squared=0.67^2
  lambda2.squared=0.67^2
  .RNG.seed<-c(1,2,3)[chain]
  .RNG.name<-c("base::Super-Duper","base::Wichmann-Hill","base::Super-Duper")
  return(list(tau=tau,preslope=preslope,thresh=thresh,
            alpha=alpha,beta=beta,alphatau=alphatau,betatau=betatau,
            lambda1.squared=lambda1.squared,lambda2.squared=lambda2.squared,
            .RNG.seed=.RNG.seed,.RNG.name=.RNG.name))
}

#Bayesian Horseshoe Model
MOD2<-"
model{
for(i in 11:N){
y[i]~dnorm(mu[i],tau)
w[i]<-1.0/(1.0+exp(-(preslope/sdy)*(y[i-2]-thresh)))
mu[i]<-(inprod(alpha[],X[i,]))*(1.0-w[i])+(inprod(beta[],X[i,]))*(w[i])
}

tau~dgamma(.001,.001)
preslope~dlnorm(3,1)
thresh~dunif(min.thresh,max.thresh)

global.squared<-global^2

for(k in 1:11){
local1.squared[k]<-(local1[k])^2
local2.squared[k]<-(local2[k])^2
alpha[k]~dnorm(0,tau/(global.squared*local1.squared[k]))
beta[k]~dnorm(0,tau/(global.squared*local2.squared[k]))
}

for(k in 1:11){
local1[k]~dt(0,1,1)T(0,)
local2[k]~dt(0,1,1)T(0,)
}

global~dt(0,1,1)T(0,)

#modules# runjags
#monitor# tau,preslope,thresh,alpha,beta
}
"

#Simulation Setting Starting value for 
initsfunction2<-function(chain){
  tau<-c(0.01,20,100)[chain]
  preslope<-c(5,10,15)[chain]
  thresh<-10.8
  set.seed(chain)
  alpha<-rnorm(11,0,1)
  global<-0.5
  local1<-rep(0.5,11)
  local2<-rep(0.5,11)
  set.seed(chain+1)
  beta<-rnorm(11,0,1)
  .RNG.seed<-c(1,2,3)[chain]
  .RNG.name<-c("base::Super-Duper","base::Wichmann-Hill","base::Super-Duper")
  return(list(tau=tau,preslope=preslope,thresh=thresh,
            alpha=alpha,beta=beta,local1=local1,local2=local2,global=global,
            .RNG.seed=.RNG.seed,.RNG.name=.RNG.name))
}

#Loop Through the Two Different Models
MOD<-list(MOD1,MOD2)
INITS<-list(initsfunction1,initsfunction2)
#If you have more than 6 cores available, change the number of clusters to 2
cl2<-makeCluster(1) 
registerDoParallel(cl2)

sunspot.out=foreach(v=1:2,.packages=c("runjags","parallel")) %dopar%{
  cl<-makeCluster(3)
  model2<-run.jags(MOD[[v]],data=datafunction(),n.chains=3,inits=INITS[[v]],
                   mutate=list(prec2sd, 'tau'),adapt=10000,
                   burnin=40000,sample=1000,thin=10,
                   method="rjparallel",method.options=list(cl=cl))
  max.psrf=max(summary(model2)[,"psrf"])
  min.ess=min(summary(model2)[,"SSeff"])
  i=1
  
  while((max.psrf>1.05|min.ess<150)&i<20){
    model2<-extend.jags(model2,adapt=1000,burnin=0,sample=1000*i,
                        silent.jags=T)
    max.psrf=max(summary(model2)[,"psrf"])
    min.ess=min(summary(model2)[,"SSeff"])
    i=i+1
    print(max.psrf)
    
  }  
  stopCluster(cl)
  out1=list(model2=model2,max.psrf=max.psrf,min.ess=min.ess,
            time=round(as.numeric(model2$timetaken)/60,1))
}
stopCluster(cl2)

###########################################################
#Plot Data From Training and Test Set (Raw and Transformed)
###########################################################
All.Data=c(Train,Test)
All.Data.transform=c(Train.transform,Test.transform)

png(file="AnnualSunspot.png",height=600,width=850)
par(mfrow=c(2,1))
plot(1700:2006,All.Data,type="l",ylab="",xlab="Year",
     main="Annual Sunspot Number")
points(1980:2006,Test,col="red",type="l")
plot(1700:2006,All.Data.transform,type="l",ylab="",xlab="Year",
     main="Square Root Transformed Annual Sunspot Number")
points(1980:2006,Test.transform,col="red",type="l")
dev.off()

###########################################################
#Get Posterior Estimates from Both Models
###########################################################
SUNSPOT.ESTIMATES<-matrix(NA,ncol=2,nrow=26)
for(k in 1:2){
  #For Bayesian Lasso
  if (k==1){
    SUNSPOT.ESTIMATES[1:25,k]=
            summary(sunspot.out[[k]]$model2)[c(4:26,2,3),"Median"] 
  }
  #For Bayesian Horseshoe
  if( k==2){
    SUNSPOT.ESTIMATES[1:25,k]=
            summary(sunspot.out[[k]]$model2)[c(4:26,2,3),"Mean"] 
  }
  SUNSPOT.ESTIMATES[26,k]=sunspot.out[[k]]$model2$sample
}




#################################################################################
#Functions Required For Recursive Forecasts 
#Using a Rolling Window Without Reestimation
#Using the Bootstrap Method for Nonlinear Model Forecasting
#################################################################################

#Function Specific For Obtaining a One Step Ahed Forecast
OneStep.func<-function(params,data,time,s=sd(Train.transform)){
  data2=c(1,data[(time-1):(time-10)])
  pred=(data2%*%params[1:11])*(1-
    (1/(1+exp(-(params[24]/s)*(data2[3]-params[25])))))+ 
    (data2%*%params[12:22])*(1/(1+exp(-(params[24]/s)*(data2[3]-params[25]))))
  return(pred)
}

#Function That Loops Through the Data Using OneStep.func for each time
#Train.Data is used to obtain residuals for 
#Bootstrapped Sampling Errors for Forecasts
MultiStep.func<-function(params,train.data,test.data,
                         s=sd(Train.transform),n.ahead){ 
  #n.ahead specifies how many time periods you would like to forecast ahead
  full.data=c(train.data,test.data,rep(NA,n.ahead))
  n.used=length(c(train.data,test.data))
  n.full=length(full.data)
  
  X.left=matrix(NA,nrow=length(train.data),ncol=10)
  for (j in 1:10){
    X.left[,j]=lag.func(train.data,k=j)
  }
  X.left=cbind(1,X.left)*(1-(1/(1+exp(-(params[24]/s)*
          (lag.func(train.data,k=2)-params[25])))))
  X.right=matrix(NA,nrow=length(train.data),ncol=10)
  for (j in 1:10){
    X.right[,j]=lag.func(train.data,k=j)
  }
  X.right=cbind(1,X.right)*(1/(1+exp(-(params[24]/s)*
            (lag.func(train.data,k=2)-params[25]))))
  X=cbind(X.left,X.right)
  predict=X%*%c(params[1:22])
  resid=train.data-predict
  resid.sample=sample(na.omit(resid),size=n.ahead,replace=T)
  
  #Forecast for long horizons using previous forecast plus random  noise
  for(i in (n.used+1):(n.used+n.ahead)){
    full.data[i]=OneStep.func(params=params,data=full.data,time=i)
    #Add randomly selected value from resampling of errors
    #from Training Data (Bootstrap Forecasts)
    full.data[i]=full.data[i]+resid.sample[i-n.used] 
  }
  return(full.data[(n.used+1):n.full])
}

#Function Used for Bootstrapping to Obtain the Pseudo-Distribution
#Of Predictions for a Specific Horizon And Can be Modified
#to Also Output Specific Forecast Quantiles
Forecast.func<-function(boot,params,train.data,test.data,
                        s=sd(Train.transform),n.ahead){
  #boot=Number of Forecasts for Each Time Period
  boot.reps=replicate(boot,MultiStep.func(params=params,
                                          train.data=train.data,
                                          test.data=test.data,
                                          n.ahead=n.ahead)) 
  #Point forecast is the mean across all Bootstrap Sampled Forecasts
  
  #Used if Forecasting One Step Ahead
  if(is.null(dim(boot.reps))) forecast=mean(boot.reps,na.rm=T)
  #Used if Forecasting More than One step Ahead
  if(!is.null(dim(boot.reps))) forecast=rowMeans(boot.reps,na.rm=T) 
  return(forecast)
}

##################################################################################
#Obtaining Forecasts Using Bayesian Lasso 
#for Horizons 1 to 5 on Transformed Test Data
#Calculating RMSFE for all Horizons by Comparing Truth to Forecasts
##################################################################################
BLASSO.FORECASTS.12345=matrix(NA,nrow=length(Test.transform),ncol=5)
for(j in 1:5){
  for(k in j:(length(Test.transform))){
    if((j-k)==0){
      BLASSO.FORECASTS.12345[k,j]=Forecast.func(boot=500,
                                          params=SUNSPOT.ESTIMATES[,1],
                                          train.data=Train.transform,
                                          test.data=NULL,
                                          n.ahead=j)[j]
    }else{
      BLASSO.FORECASTS.12345[k,j]=Forecast.func(boot=500,
                                          params=SUNSPOT.ESTIMATES[,1],
                                          train.data=Train.transform,
                                          test.data=Test.transform[1:(k-j)],
                                          n.ahead=j)[j]
    }
  }  
} 

RMSFE1=rep(NA,5)
for(k in 1:5){
  #Other Metrics for Forecast Evaluation Can be Replaced Here
  RMSFE1[k]=sqrt(mean((Test.transform-BLASSO.FORECASTS.12345[,k])^2,na.rm=T)) 
}

##################################################################################
#Obtaining Forecasts Using Bayesian Horseshoe 
#for Horizons 1 to 5 on Transformed Test Data
#Calculating RMSFE for all Horizons by Comparing Truth to Forecasts
##################################################################################
BHS.FORECASTS.12345=matrix(NA,nrow=length(Test.transform),ncol=5)
for(j in 1:5){
  for(k in j:(length(Test.transform))){
    if((j-k)==0){
      BHS.FORECASTS.12345[k,j]=Forecast.func(boot=500,
                                             params=SUNSPOT.ESTIMATES[,2],
                                             train.data=Train.transform,
                                             test.data=NULL,
                                             n.ahead=j)[j]
    }else{
      BHS.FORECASTS.12345[k,j]=Forecast.func(boot=500,
                                             params=SUNSPOT.ESTIMATES[,2],
                                             train.data=Train.transform,
                                             test.data=Test.transform[1:(k-j)],
                                             n.ahead=j)[j]
    }
  }  
} 

RMSFE2=rep(NA,5)
for(k in 1:5){
  RMSFE2[k]=sqrt(mean((Test.transform-BHS.FORECASTS.12345[,k])^2,na.rm=T))
}

