library(mitools)
library(survey)

system.time({
ccr<-replicate(1000,{
    df<-data.frame(y=rep(1:0,c(500,10000-500)))
    df$x<-with(df, rnorm(10000,y*3-1.5))

### Case-control sampling
insample<-c(which(df$y==1),sample(which(df$y==0),sum(df$y)))
sdf<-df[insample,]
fdf<-df
fdf$x[-insample]<-NA
Nmiss<-nrow(df)-length(insample)
w<-sum(df$y==0)/sum(sdf$y==0)
sdf$wt<-with(sdf,ifelse(y==1,1,w))

#fits
model0<-glm(y~x,family=binomial,data=df)
modelm<-glm(y~x,family=binomial,data=sdf)
models<-glm(y~x,family=quasibinomial,data=sdf,weight=wt)


### Resampling imputation
p3<-vector("list",10)
for(i in 1:10){
  p3[[i]]<-fdf
  p3[[i]]$x[-insample]<-with(sdf,sample(x[y==0],Nmiss,replace=TRUE))
}
p3list<-imputationList(p3)

# Parametric imputation
p4<-vector("list",10)
    m0<-with(sdf, mean(x-y*coef(modelm)[2]))
    s0<-with(sdf, sd(x-y*coef(modelm)[2]))
    
for(i in 1:10){
  p4[[i]]<-fdf
  p4[[i]]$x[-insample]<-rnorm(Nmiss,m0,s0)
}
p4list<-imputationList(p4)




## fit
modelmi1<-with(p3list, glm(y~x,family=binomial))
modelmi2<-with(p4list, glm(y~x,family=binomial))

beta<-c(coef((modelm)),
        coef((models)),
        coef(MIcombine(modelmi1)),
        coef(MIcombine(modelmi2))
       )-coef(model0)

beta
}
)
})
