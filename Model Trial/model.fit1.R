data(gutten)
tree1.1<-gutten[gutten$tree.ID=="1.1",]
plot(tree1.1$age.bh,tree1.1$dbh.cm,type='p')
model.fit1<-nls(dbh.cm~B1*(1-exp(-log(2)/B2*age.bh)),
                start = list(B1=29,B2=10),
                data = tree1.1)
summary(model.fit1)

#function
dbh.pred<-function(age.bh,B1,B2){
  dbh.p<-B1*(1-exp(-log(2)/B2*age.bh))
  return(dbh.p)}

model.fit2<-nls(dbh.cm~dbh.pred(age.bh,B1,B2),
                start=list(B1=29,B2=10),
                data=tree1.1)
summary(model.fit2)

tval<-qt(0.975,summary(model.fit1)$df[2])
coef(summary(model.fit1))[,1:2]%*%matrix(c(1,-tval,1,tval),nrow=2)

library(nlstools)
