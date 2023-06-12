setwd("/Users/mingtzu/Rcode")
library("nlstools")
library("nlme")
library(lattice)
library(lmtest)

gutten2<-read.csv("gutten2.csv")
gutten2$BA<-(gutten2$dbh.cm)^2*0.00007854
for (i in 1:1200){
  gutten2$num[i]=paste(gutten2$location[i],gutten2$tree[i],sep=".")
  sumBA<-0
  for (j in 1:1200){
    if(gutten2$age.base[j]==gutten2$age.base[i] & gutten2$location[j]==gutten2$location[i]){
      if (gutten2$BA[j]>gutten2$BA[i]){
        sumBA<-sumBA+gutten2$BA[j]
      }
    }
  }
  gutten2$BAL[i]<-sumBA
}
#View(gutten2)
#Primitive model of Ber NLME
dbh.pred1<-function(age.base,B1,B2){
  dbh.p<-(B1)*(1-exp(-(B2)*age.base))^3
  return(dbh.p)}

nls.fit.Ber1<-nls(dbh.cm~dbh.pred1(age.base,B1,B2),
                  start=list(B1=27.360790,B2=0.054481),
                  data=gutten2)
summary(nls.fit.Ber1)
plot(nlsResiduals(nls.fit.Ber1),which=1)
gutten2.d.Ber1<-groupedData(dbh.cm~age.base|num,data=gutten2)
gutten2.nlsList.Ber1<-nlsList(dbh.cm~dbh.pred1(age.base,B1,B2),
                              start=list(B1=27,B2=0.05),
                              data=gutten2.d.Ber1,subset = gutten2.d.Ber1$location%in%c(1,3,4,5,7))
summary(gutten2.nlsList.Ber1)
gutten2.nlme.Ber1<-nlme(gutten2.nlsList.Ber1)
summary(gutten2.nlme.Ber1)
gutten2.nlme.Ber1<-update(gutten2.nlme.Ber1,correlation= corARMA (p=1,q=0),control=list( maxIter =200))
plot(ACF(gutten2.nlme.Ber1,resType='n',maxLag=10),alpha=0.01)

#full model of Ber NLME
dbh.pred2<-function(age.base,si,BAL,B10,B11,B12,B20,B21,B22){
  dbh.p<-(B10+B11*si+B12*BAL)*(1-exp(-(B20+B21*si+B22*BAL)*age.base))^3
  return(dbh.p)}

nls.fit.Ber2<-nls(dbh.cm~dbh.pred2(age.base,si,BAL,B10,B11,B12,B20,B21,B22),
                  start=list(B10=27.360790,B11=0,B12=0,B20=0.054481,B21=0,B22=0),
                  data=gutten2)
summary(nls.fit.Ber2)
# gutten2.d.Ber2<-groupedData(dbh.cm~age.base|num,data=gutten2)
# gutten2.nlme.Ber2<-nlme(dbh.cm~dbh.pred2(age.base,si,BAL,B10,B11,B12,B20,B21,B22),
#                         start=c(B10=27.360790,B11=0,B12=0,B20=0.054481,B21=0,B22=0),
#                         data=gutten2.d.Ber2,
#                         fixed=B10+B11+B12+B20+B21+B22~1,
#                         subset=gutten2.d.Ber2$location%in%c(1,3,4,5,7))
# summary(gutten2.nlme.Ber2)
# 
# gutten2.nlme.Ber2<-update(gutten2.nlme.Ber2,correlation= corARMA (p=1,q=0),control=list( maxIter =200))
# plot(ACF(gutten2.nlme.Ber2,resType='n',maxLag=10),alpha=0.01)
#多行註解：command+shift+c



gutten2.d.Ber2<-groupedData(dbh.cm~age.base|location/num,data=gutten2)
gutten2.nlme.Ber2<-nlme(dbh.cm~dbh.pred2(age.base,si,BAL,B10,B11,B12,B20,B21,B22),
                        start=c(B10=27.360790,B11=0,B12=0,B20=0.054481,B21=0,B22=0),
                        data=gutten2.d.Ber2,
                        fixed=B10+B11+B12+B20+B21+B22~1,
                        random=B10+B20~1,subset=gutten2.d.Ber2$location%in%c(1,3,4,5,7))
summary(gutten2.nlme.Ber2)
gutten2.nlme.Ber2<-update(gutten2.nlme.Ber2,correlation= corARMA (p=1,q=0),control=list( maxIter =200))
plot(ACF(gutten2.nlme.Ber2,resType='n',maxLag=10),alpha=0.01)



#Primitive model of Gom NLME
dbh.pred_1<-function(age.base,B1,B2,B3){
  dbh.p<-B1*(exp(-B2*exp(-B3*age.base)))
  return(dbh.p)}
nls.fit.Gom1<-nls(dbh.cm~dbh.pred_1(age.base,B1,B2,B3),
                  start=list(B1=28.715611,B2=2.343944,B3=0.039787),
                  data=gutten2)
summary(nls.fit.Gom1)
plot(nlsResiduals(nls.fit.Gom1),which=1)
gutten2.d.Gom1<-groupedData(dbh.cm~age.base|num,data=gutten2)
gutten2.nlsList.Gom1<-nlsList(dbh.cm~dbh.pred_1(age.base,B1,B2,B3),
                              start=list(B1=28.715611,B2=2.343944,B3=0.039787),
                              data=gutten2.d.Gom1,subset = gutten2.d.Gom1$location%in%c(1,3,4,5,7))
summary(gutten2.nlsList.Gom1)
gutten2.nlme.Gom1<-nlme(gutten2.nlsList.Gom1)
summary(gutten2.nlme.Gom1)
gutten2.nlme.Gom1<-update(gutten2.nlme.Gom1,correlation= corARMA (p=0,q=1),control=list( maxIter =200))
plot(ACF(gutten2.nlme.Gom1,resType='n',maxLag=10),alpha=0.01)

#full model of Gom NLME
dbh.pred_2<-function(age.base,si,BAL,B10,B11,B12,B20,B21,B30,B31){
  dbh.p<-(B10+B11*si+B12*BAL)*(exp(-(B20+B21*si)*exp(-(B30+B31*si)*age.base)))
  return(dbh.p)}
nls.fit.Gom2<-nls(dbh.cm~dbh.pred_2(age.base,si,BAL,B10,B11,B12,B20,B21,B30,B31),
                  start=list(B10=28.715611,B11=0,B12=0,B20=2.343944,B21=0,B30=0.039787,B31=0),
                  data=gutten2)
summary(nls.fit.Gom2)
# gutten2.d.Gom2<-groupedData(dbh.cm~age.base|num,data=gutten2)
# gutten2.nlme.Gom2<-nlme(dbh.cm~dbh.pred_2(age.base,si,BAL,B10,B11,B12,B20,B21,B22,B30,B31),
#                         start=c(B10=28.715611,B11=0,B12=0,B20=2.343944,B21=0,B22=0,B30=0.039787,B31=0),
#                         data=gutten2.d.Gom2,
#                         fixed=B10+B11+B12+B20+B21+B22+B30+B31~1,
#                         subset=gutten2.d.Gom2$location%in%c(1,3,4,5,7))
# summary(gutten2.nlme.Gom2)
# 
# gutten2.nlme.Gom2<-update(gutten2.nlme.Gom2,correlation= corARMA (p=1,q=0),control=list( maxIter =200))
# plot(ACF(gutten2.nlme.Gom2,resType='n',maxLag=10),alpha=0.01)

gutten2.d.Gom2<-groupedData(dbh.cm~age.base|location/num,data=gutten2)
gutten2.nlme.Gom2<-nlme(dbh.cm~dbh.pred_2(age.base,si,BAL,B10,B11,B12,B20,B21,B30,B31),
                        start=c(B10=28.715611,B11=0,B12=0,B20=2.343944,B21=0,B30=0.039787,B31=0),
                        data=gutten2.d.Gom2,
                        fixed=B10+B11+B12+B20+B21+B30+B31~1,
                        random=B10+B20+B30~1,subset=gutten2.d.Gom2$location%in%c(1,3,4,5,7))
summary(gutten2.nlme.Gom2)
gutten2.nlme.Gom2<-update(gutten2.nlme.Gom2,correlation= corARMA (p=1,q=0),control=list( maxIter =200))
plot(ACF(gutten2.nlme.Gom2,resType='n',maxLag=10),alpha=0.01)
lrtest(gutten2.nlme.Ber2,gutten2.nlme.Gom2)
AIC(gutten2.nlme.Ber2)
AIC(gutten2.nlme.Gom2)
