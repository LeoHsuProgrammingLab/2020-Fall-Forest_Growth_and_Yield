plot(fitted(model.fit1),tree1.1$dbh.cm,xlab="
Fitted Values",ylab="Observed Values")
abline(0,1,col="red")

plot(nlsResiduals(model.fit1),which=1)#Raw residuals
plot(nlsResiduals(model.fit1),which=2)#Standard residuals沒有單位，只有比較的值
plot(nlsResiduals(model.fit1),which=6)#QQplot，如果是normally distribution，點都會在線上
vcov(model.fit1)#co-e=cov/(sigma(a)*sigma(b))
vcov(model.fit1)[1,2]/prod(summary(model.fit1)$coefficients[1,2])#co-efficients

