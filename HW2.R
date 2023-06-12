library(FAwR)
data(gutten)
tree1.1<-gutten[gutten$tree.ID=="1.1",]
plot(tree1.1$age.bh,tree1.1$dbh.cm,type='p')

#Hossfeld
dbh.pred<-function(age.bh,B1,B2,B3){
  dbh.p<-age.bh^B3/((age.bh)^B3/B1+B2)
  return(dbh.p)}
curve(dbh.pred(x,33.81748,4.35425,1.39638),from=0,to=200)
points(tree1.1$age.bh,tree1.1$dbh.cm)

model.fit.Hos<-nls(dbh.cm~dbh.pred(age.bh,B1,B2,B3),
                start=list(B1=33.81748,B2=4.35425,B3=1.39638),
                data=tree1.1)
summary(model.fit.Hos)
vcov(model.fit.Hos)[1,2]/prod(summary(model.fit.Hos)$coefficients[,2])
dbh.pred(200,33.81748,4.35425, 1.39638)

#Bertalanffy
dbh.pred<-function(age.bh,B1,B2){
  dbh.p<-B1*(1-exp(-B2*age.bh))^3
  return(dbh.p)}
curve(dbh.pred(x,27.360790,0.054481),from=0,to=200)
points(tree1.1$age.bh,tree1.1$dbh.cm)

model.fit.Ber<-nls(dbh.cm~dbh.pred(age.bh,B1,B2),
                start=list(B1=27.360790,B2=0.054481),
                data=tree1.1)
summary(model.fit.Ber)
dbh.pred(200,27.360790,0.054481)

#Chapman
dbh.pred<-function(age.bh,B1,B2,B3){
  dbh.p<-B1*(1-exp(-B2*age.bh))^B3
  return(dbh.p)}
curve(dbh.pred(x,29.879787,0.026848,1.205632),from=0,to=200)
points(tree1.1$age.bh,tree1.1$dbh.cm)

model.fit.Cha<-nls(dbh.cm~dbh.pred(age.bh,B1,B2,B3),
                start=list(B1=30,B2=0.04,B3=3),
                data=tree1.1)
summary(model.fit.Cha)
dbh.pred(200,29.879787,0.026848,1.205632)

#Gompertz
dbh.pred<-function(age.bh,B1,B2,B3){
  dbh.p<-B1*(exp(-B2*exp(-B3*age.bh)))
  return(dbh.p)}
curve(dbh.pred(x,28.715611,2.343944,0.039787),from=0,to=200)
points(tree1.1$age.bh,tree1.1$dbh.cm)

model.fit.Gom<-nls(dbh.cm~dbh.pred(age.bh,B1,B2,B3),
                start=list(B1=30,B2=30,B3=0.12),
                data=tree1.1)
summary(model.fit.Gom)
dbh.pred(200,28.715611,2.343944, 0.039787)

#Weibull
dbh.pred<-function(age.bh,B1,B2,B3){
  dbh.p<-B1*(1-exp(-B2*age.bh^B3))
  return(dbh.p)}
curve(dbh.pred(x,29.716360,0.015413,1.111638),from=0,to=200)
points(tree1.1$age.bh,tree1.1$dbh.cm)

model.fit.Wei<-nls(dbh.cm~dbh.pred(age.bh,B1,B2,B3),
                start=list(B1=30,B2=0.025,B3=1),
                data=tree1.1)
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
B1_B2<-vcov(model.fit.Hos)[1,2]/prod(summary(model.fit.Hos)$coefficients[1,2],summary(model.fit.Hos)$coefficients[2,2])
B2_B3<-vcov(model.fit.Hos)[2,3]/prod(summary(model.fit.Hos)$coefficients[2,2],summary(model.fit.Hos)$coefficients[3,2])
B1_B3<-vcov(model.fit.Hos)[1,3]/prod(summary(model.fit.Hos)$coefficients[1,2],summary(model.fit.Hos)$coefficients[3,2])

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
qt(0.975,summary(model.fit.Hos)$df[2])*0.25547
qt(0.975,summary(model.fit.Hos)$df[2])*0.26870
qt(0.975,summary(model.fit.Hos)$df[2])*0.01942
summary(model.fit.Hos)
qt(0.975,summary(model.fit.Ber)$df[2])*0.681044
qt(0.975,summary(model.fit.Ber)$df[2])*0.003401
summary(model.fit.Ber)
qt(0.975,summary(model.fit.Cha)$df[2])*0.293698
qt(0.975,summary(model.fit.Cha)$df[2])*0.001308
qt(0.975,summary(model.fit.Cha)$df[2])*0.052280
summary(model.fit.Cha)
qt(0.975,summary(model.fit.Gom)$df[2])*0.444456
qt(0.975,summary(model.fit.Gom)$df[2])*0.163924
qt(0.975,summary(model.fit.Gom)$df[2])*0.002657

  
summary(model.fit.Gom)

qt(0.975,summary(model.fit.Wei)$df[2])*0.372122
qt(0.975,summary(model.fit.Wei)$df[2])*0.001624
qt(0.975,summary(model.fit.Wei)$df[2])*0.032728

summary(model.fit.Wei)


#Project each of the five fitted models to 200 years

