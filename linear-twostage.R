library(missreg3)
library(survey)
library(splines)
library(mitools)


estfun.glm<-function (model) 
{
    xmat <- model.matrix(model)
    residuals(model, "working") * model$weights * xmat
}


epsilon<-0.05
one.sim<-function(epsilon){
df<-data.frame(x=rnorm(4000))
df$z<-df$x+rnorm(4000)
df$zstrat<-cut(df$z, c(-Inf, qnorm(.05,s=sqrt(2)),qnorm(0.95,s=sqrt(2)),Inf))
df$zmidx<-ifelse(as.numeric(df$zstrat)==2,df$x,0)
df$y<-df$x+rnorm(4000)-epsilon*df$zmidx
df$id<-1:nrow(df)
df$fullx<-df$x

df$insample<-df$id %in% c(which(df$zstrat=="(-Inf,-2.33]"), which(df$zstrat=="(2.33, Inf]"),sample(which(df$zstrat=="(-2.33,2.33]"),200))
df$x[!df$insample]<-NA

impmodel2<-glm(x~y+z,data=subset(df,insample))


## single imputation
df$predx2<-predict(impmodel2,newdata=df)

calmodel2<-glm(y~predx2,data=df)

eff<-as.data.frame(estfun.glm(calmodel2))
names(eff)<-c("eff1","eff2")

df2<-cbind(df,eff)

## uncalibrated
des<-twophase(id=list(~id,~id), strata=list(NULL,~zstrat), data=df2, subset=~insample)

## calibrated
cdes2<-calibrate(des,formula=~(eff1+eff2)*zstrat,calfun="raking")

## MI
p<-vector("list",10)
    
for(i in 1:10){
    p[[i]]<-df
    pred<-predict(impmodel2, newdata=df,se.fit=TRUE)
    p[[i]]$predx<- rnorm(nrow(df), pred$fit, sqrt(pred$se.fit^2+pred$residual.scale^2))
}
plist<-imputationList(p)
modelmi<-with(plist, glm(y~predx))
effs<-as.data.frame(Reduce("+",lapply(modelmi, estfun.glm))/10)
names(effs)<-c("mieff1","mieff2")
df3<-cbind(df,effs)

des<-twophase(id=list(~id,~id), strata=list(NULL,~zstrat), data=df3, subset=~insample)
## calibrated
cdesmi<-calibrate(des,formula=~(mieff1+mieff2)*zstrat,calfun="raking")


## mle
df$obstype.name<-ifelse(df$insample,"retro","strata")
yCuts<-c(-3,-2.33,-2,-1,0,1,2,2.33,3)
mz<-locsc2stg(y~x,~1,xstrata="zstrat",data=df,xs.includes=FALSE,method="ycutmeth",obstype.name="obstype.name",yCuts=yCuts,errdistn="normal", print.progress=0)
mz2<-locsc2stg(y~x,~1,xstrata="zstrat",data=df,xs.includes=FALSE,method="direct",obstype.name="obstype.name",start=mz$coefficients,errdistn="normal", print.progress=0)


cat(".")
flush.console()
c(
uncal=coef(svyglm(y~x,design=des))[2],
impcal=coef(svyglm(y~x,design=cdes2))[2],
mical=coef(svyglm(y~x,design=cdesmi))[2],
    mi=coef(MIcombine(modelmi))[2],
mle=summary(mz2)$coef.table2[2,1],
census=coef(glm(y~fullx,data=df))[2]
)}

epsilons<-c(0,0.025,0.05,0.06,0.07,0.08,0.1,0.15)
rr2<-lapply(epsilons, function(eps) {print(eps);replicate(1000, one.sim(eps))})
save(rr2,epsilons,file="linear-twostage.rda")
