library(FAwR)
library(nlstools)
library(lmtest)
data(gutten)
plot(gutten$age.bh,gutten$dbh.cm,type='p')

#Hossfeld
dbh.pred<-function(age.bh,B1,B2,B3){
  dbh.p<-age.bh^B3/((age.bh)^B3/B1+B2)
  return(dbh.p)}
curve(dbh.pred(x,33.81748,4.35425,1.39638),from=0,to=200,ylim=c(0,max(gutten$dbh.cm)))
points(gutten$age.bh,gutten$dbh.cm)

model.fit.Hos<-nls(dbh.cm~dbh.pred(age.bh,B1,B2,B3),
                   start=list(B1=33.81748,B2=4.35425,B3=1.39638),
                   data=gutten)
summary(model.fit.Hos)
vcov(model.fit.Hos)[1,2]/prod(summary(model.fit.Hos)$coefficients[,2])
dbh.pred(200,33.81748,4.35425, 1.39638)

#Bertalanffy
dbh.pred<-function(age.bh,B1,B2){
  dbh.p<-B1*(1-exp(-B2*age.bh))^3
  return(dbh.p)}
curve(dbh.pred(x,27.360790,0.054481),from=0,to=200,ylim=c(0,max(gutten$dbh.cm)))
points(gutten$age.bh,gutten$dbh.cm)

model.fit.Ber<-nls(dbh.cm~dbh.pred(age.bh,B1,B2),
                   start=list(B1=27.360790,B2=0.054481),
                   data=gutten)
summary(model.fit.Ber)
dbh.pred(200,27.360790,0.054481)

#Chapman
dbh.pred<-function(age.bh,B1,B2,B3){
  dbh.p<-B1*(1-exp(-B2*age.bh))^B3
  return(dbh.p)}
curve(dbh.pred(x,29.879787,0.026848,1.205632),from=0,to=200,ylim=c(0,max(gutten$dbh.cm)))
points(gutten$age.bh,gutten$dbh.cm)

model.fit.Cha<-nls(dbh.cm~dbh.pred(age.bh,B1,B2,B3),
                   start=list(B1=30,B2=0.04,B3=3),
                   data=gutten)
summary(model.fit.Cha)
dbh.pred(200,29.879787,0.026848,1.205632)

#Gompertz
dbh.pred<-function(age.bh,B1,B2,B3){
  dbh.p<-B1*(exp(-B2*exp(-B3*age.bh)))
  return(dbh.p)}
curve(dbh.pred(x,28.715611,2.343944,0.039787),from=0,to=200,ylim=c(0,max(gutten$dbh.cm)))
points(gutten$age.bh,gutten$dbh.cm)

model.fit.Gom<-nls(dbh.cm~dbh.pred(age.bh,B1,B2,B3),
                   start=list(B1=40,B2=30,B3=0.12),
                   data=gutten)
summary(model.fit.Gom)
dbh.pred(200,28.715611,2.343944, 0.039787)

#Weibull
dbh.pred<-function(age.bh,B1,B2,B3){
  dbh.p<-B1*(1-exp(-B2*age.bh^B3))
  return(dbh.p)}
curve(dbh.pred(x,29.716360,0.015413,1.111638),from=0,to=200,ylim=c(0,max(gutten$dbh.cm)))
points(gutten$age.bh,gutten$dbh.cm)

model.fit.Wei<-nls(dbh.cm~dbh.pred(age.bh,B1,B2,B3),
                   start=list(B1=30,B2=0.025,B3=1),
                   data=gutten)
summary(model.fit.Wei)
dbh.pred(200,29.716360,0.015413, 1.111638)

#Compare the residuals
plot(nlsResiduals(model.fit.Hos),which=1)
plot(nlsResiduals(model.fit.Ber),which=1)
plot(nlsResiduals(model.fit.Cha),which=1)
plot(nlsResiduals(model.fit.Gom),which=1)
plot(nlsResiduals(model.fit.Wei),which=1)

#Compare goodness of fit
plot(nlsResiduals(model.fit.Hos),which=2)
plot(nlsResiduals(model.fit.Ber),which=2)
plot(nlsResiduals(model.fit.Cha),which=2)
plot(nlsResiduals(model.fit.Gom),which=2)
plot(nlsResiduals(model.fit.Wei),which=2)

#Compare correlation between fitted parameters of each of the five model forms
vcov(model.fit.Hos)[1,2]/prod(summary(model.fit.Hos)$coefficients[1,2],summary(model.fit.Hos)$coefficients[2,2])
vcov(model.fit.Hos)[2,3]/prod(summary(model.fit.Hos)$coefficients[2,2],summary(model.fit.Hos)$coefficients[3,2])
vcov(model.fit.Hos)[1,3]/prod(summary(model.fit.Hos)$coefficients[1,2],summary(model.fit.Hos)$coefficients[3,2])

vcov(model.fit.Ber)[1,2]/prod(summary(model.fit.Ber)$coefficients[1,2],summary(model.fit.Ber)$coefficients[2,2])

vcov(model.fit.Cha)[1,2]/prod(summary(model.fit.Cha)$coefficients[1,2],summary(model.fit.Cha)$coefficients[2,2])
vcov(model.fit.Cha)[2,3]/prod(summary(model.fit.Cha)$coefficients[2,2],summary(model.fit.Cha)$coefficients[3,2])
vcov(model.fit.Cha)[1,3]/prod(summary(model.fit.Cha)$coefficients[1,2],summary(model.fit.Cha)$coefficients[3,2])

vcov(model.fit.Gom)[1,2]/prod(summary(model.fit.Gom)$coefficients[1,2],summary(model.fit.Gom)$coefficients[2,2])
vcov(model.fit.Gom)[2,3]/prod(summary(model.fit.Gom)$coefficients[2,2],summary(model.fit.Gom)$coefficients[3,2])
vcov(model.fit.Gom)[1,3]/prod(summary(model.fit.Gom)$coefficients[1,2],summary(model.fit.Gom)$coefficients[3,2])

vcov(model.fit.Wei)[1,2]/prod(summary(model.fit.Wei)$coefficients[1,2],summary(model.fit.Wei)$coefficients[2,2])
vcov(model.fit.Wei)[2,3]/prod(summary(model.fit.Wei)$coefficients[2,2],summary(model.fit.Wei)$coefficients[3,2])
vcov(model.fit.Wei)[1,3]/prod(summary(model.fit.Wei)$coefficients[1,2],summary(model.fit.Wei)$coefficients[3,2])


#Compute 95% confidence intervals
qt(0.975,summary(model.fit.Hos)$df[2])*summary(model.fit.Hos)$coefficients[1,2]
qt(0.975,summary(model.fit.Hos)$df[2])*summary(model.fit.Hos)$coefficients[2,2]
qt(0.975,summary(model.fit.Hos)$df[2])*summary(model.fit.Hos)$coefficients[3,2]
summary(model.fit.Hos)

qt(0.975,summary(model.fit.Ber)$df[2])*summary(model.fit.Ber)$coefficients[1,2]
qt(0.975,summary(model.fit.Ber)$df[2])*summary(model.fit.Ber)$coefficients[2,2]
summary(model.fit.Ber)

qt(0.975,summary(model.fit.Cha)$df[2])*summary(model.fit.Cha)$coefficients[1,2]
qt(0.975,summary(model.fit.Cha)$df[2])*summary(model.fit.Cha)$coefficients[2,2]
qt(0.975,summary(model.fit.Cha)$df[2])*summary(model.fit.Cha)$coefficients[3,2]
summary(model.fit.Cha)

qt(0.975,summary(model.fit.Gom)$df[2])*summary(model.fit.Gom)$coefficients[1,2]
qt(0.975,summary(model.fit.Gom)$df[2])*summary(model.fit.Gom)$coefficients[2,2]
qt(0.975,summary(model.fit.Gom)$df[2])*summary(model.fit.Gom)$coefficients[3,2]
summary(model.fit.Gom)

qt(0.975,summary(model.fit.Wei)$df[2])*summary(model.fit.Wei)$coefficients[1,2]
qt(0.975,summary(model.fit.Wei)$df[2])*summary(model.fit.Wei)$coefficients[2,2]
qt(0.975,summary(model.fit.Wei)$df[2])*summary(model.fit.Wei)$coefficients[3,2]
summary(model.fit.Wei)

dbh.pred<-function(age.bh,B1,B2,B3){
  dbh.p<-B1*(1-B2*exp(B3*age.bh))
  return(dbh.p)}
curve(dbh.pred(x,25.782333333333337,1.1129999999999998,-0.02766666666666667),from=0,to=200,ylim=c(0,max(gutten$dbh.cm)))
points(gutten$age.bh,gutten$dbh.cm)

dbh.pred<-function(age.bh,B1,B2,B3){
  dbh.p<-exp(B1*(1+B2*exp(B3*age.bh)))
  return(dbh.p)}
curve(dbh.pred(x,5.692150000000001,0.5702000000000002,-0.026500000000000006),from=0,to=200)
points(gutten$age.bh,gutten$dbh.cm)

#Project each of the five fitted models to 200 years
lrtest(model.fit.Hos,model.fit.Ber)#越高越好
lrtest(model.fit.Hos,model.fit.Cha)
lrtest(model.fit.Hos,model.fit.Gom)
lrtest(model.fit.Hos,model.fit.Wei)
lrtest(model.fit.Ber,model.fit.Cha)
lrtest(model.fit.Ber,model.fit.Gom)
lrtest(model.fit.Ber,model.fit.Wei)
lrtest(model.fit.Cha,model.fit.Gom)
lrtest(model.fit.Cha,model.fit.Wei)
lrtest(model.fit.Gom,model.fit.Wei)
AIC(model.fit.Hos)#越低越好
AIC(model.fit.Ber)
AIC(model.fit.Cha)
AIC(model.fit.Gom)
AIC(model.fit.Wei)


