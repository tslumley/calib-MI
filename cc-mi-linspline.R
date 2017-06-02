library(mitools)
library(survey)

expit<-function(eta) exp(eta)/(1+exp(eta))



makeData<-function(beta=0.07,knot=1.8){
df<-data.frame(x=rnorm(1e5))
df$mu<-expit(df$x-5+beta*pmax(df$x-knot,0))
df$y<-rbinom(nrow(df),1,df$mu)

cases<-which(df$y==1)
controls<-sample(which(df$y==0),length(cases))

df$insample<-(1:nrow(df)) %in% c(cases,controls)


df$wt<-ifelse(df$y==1,1,sum(df$y==0)/length(cases))
df$wt[!df$insample]<-NA
df
}

doCensus<-function(df){
c(b=coef(glm(y~x,data=df, family=binomial))[2])
}

doMLE<-function(df){
m1<-glm(y~x,data=df,subset=insample, family=binomial)
c(b=coef(m1)[2],s=SE(m1)[2])
}

doSurvey<-function(df){
des2<-twophase(id=list(~1,~1), strata=list(NULL,~y), data=df, method="simple",subset=~insample)
m2<-svyglm(y~x,design=des2,family=quasibinomial)
c(b=coef(m2)[2],s=SE(m2)[2])
}

doMIresample<-function(df){
Nmiss<-sum(!df$insample)
### Resampling imputation
p3<-vector("list",10)
for(i in 1:10){
  p3[[i]]<-df
  p3[[i]]$x[!df$insample]<-with(subset(df,insample),sample(x[y==0],Nmiss,replace=TRUE))
}
p3list<-imputationList(p3)
modelmi1<-MIcombine(with(p3list, glm(y~x,family=binomial)))
c(b=coef(modelmi1)[2], s=SE(modelmi1)[2])
}

doMINormal<-function(df){
    Nmiss<-sum(!df$insample)
    modelm<-lm(x~y,data=df,subset=insample)
# Parametric imputation
p4<-vector("list",10)
    m0<-with(subset(df,insample), mean(x-y*coef(modelm)[2]))
    s0<-with(subset(df,insample), sd(x-y*coef(modelm)[2]))
    
for(i in 1:10){
  p4[[i]]<-df
  p4[[i]]$x[!df$insample]<-rnorm(Nmiss,m0,s0)
}
p4list<-imputationList(p4)
modelmi2<-MIcombine(with(p4list, glm(y~x,family=binomial)))
c(b=coef(modelmi2)[2], s=SE(modelmi2)[2])
}



simOne<-function(beta=0.3){
    d1<-makeData(beta)
    c(census=doCensus(d1),
      survey=doSurvey(d1),
      mle=doMLE(d1),
      MIr=doMIresample(d1),
      MIg=doMINormal(d1)
      )
      
}


sims<-lapply(c(0,0.2,0.3,0.35,0.4,0.6), function(b) {print(b); replicate(1000, simOne(b))})
save(sims,file="cc-mi-linspline.rda")


