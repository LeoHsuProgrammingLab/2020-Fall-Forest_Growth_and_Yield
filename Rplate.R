setwd("/Users/mingtzu/Rcode")
library("nlstools")
library("nlme")
library(lattice)

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
View(gutten2)
#Hossfeld的nls fit
dbh.pred<-function(age.base,si,B10,B11,B20,B21,B30,B31){
  dbh.p<-age.base^(B30+B31*si)/((age.base^(B30+B31*si)/(B10+B11*si))+(B20+B21*si))
  return(dbh.p)}
nls.fit.Hos<-nls(dbh.cm~dbh.pred(age.base,si,B10,B11,B20,B21,B30,B31),
                 start=list(B10=33,B11=0,B20=4,B21=0,B30=1,B31=0),
                 data=gutten2)
summary(nls.fit.Hos)
plot(nlsResiduals(nls.fit.Hos),which=1)

#nlme
gutten2.d<-groupedData(dbh.cm~age.base|num,data=gutten2)
head(gutten2.d)
str(gutten2.d)

gutten2.nls<-nls(dbh.cm~dbh.pred(age.base,si,B10,B11,B20,B21,B30,B31),
                 start=list(B10=33,B11=0,B20=4,B21=0,B30=1,B31=0),
                 data=gutten2)
summary(gutten2.nls)
gutten2.nlsList<-nlsList(dbh.cm~dbh.pred(age.base,si,B10,B11,B20,B21,B30,B31),
                         start=list(B10=33,B11=0,B20=4,B21=0,B30=1,B31=0),
                         data=gutten2.d)
summary(gutten2.nlsList)
gutten2.nlme.0<-nlme(gutten2.nlsList)
summary(gutten2.nlme.0)