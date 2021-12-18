## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----echo = TRUE, eval = TRUE-------------------------------------------------
x <- rnorm(200)
fig <- hist(x, breaks=10, col='blue')

## ----include=FALSE------------------------------------------------------------
x <- c(1:500)
y <- rnorm(500,x,1)
b <- lm(y ~ x)
df <- summary(b)$coef
print(df)

## ----echo=FALSE---------------------------------------------------------------
knitr::kable(df)


## ----eval=TRUE----------------------------------------------------------------
   Yray <- function(z,sigma){y<-((2*sigma^2)*log(1/(1-z)))^(1/2);return(y)} ##cdf反函数
   Ans <- function(n,sigma){
     u <- runif(n);
     x<-Yray(u,sigma);
     hist(x, prob = TRUE, main = expression(f(x)==(x/σ^2)*exp((-x^2)/(2*σ^2))));
     y <- seq(0, 4*sigma,0.01);
     lines(y, (y/sigma^2)*exp((-y^2)/(2*sigma^2)))
     } ## n代表取的随机数个数
   Ans(10000,0.05)
   Ans(10000,0.5)
   Ans(10000,5)
   Ans(10000,50)

## ----eval=TRUE----------------------------------------------------------------
    Ans2 <- function(p){
    u<-rnorm(1000);
    for(i in 1:1000){
      if(runif(1)>p){
        u[i]<-u[i]+3}};
    return(u)
    } ##生成1000个按要求的数，接下来对p以0.01为间隔生成多个直方图
    for(i in 1:20){hist(Ans2(i/20),breaks=20,prob=TRUE,main='',xlab='')} ## p从0.01，0.02一直到0.10，下类似
    ##依据所得的结果猜测[0.35,0.65]区间内有明显的双峰效应

## ----eval=TRUE----------------------------------------------------------------
    Ans3 <- function(n,t,lambda,alpha,beta){
      X<-rep(0,n);
      for(i in 1:n){
        X[i]<-sum(rgamma(rpois(1,t*lambda),alpha,beta))};
      return(X)
      } ##生成n个这样的X(t),E(X(t))理论值lambda*t*alpha/beta,VAR(X(t))理论值lambda*t*(alpha+alpha^2)/beta.
    mv <- function(n,t,lambda,alpha,beta){
      mv<-rep(0,4);
      mv[1]<-mean(Ans3(n,t,lambda,alpha,beta));
      mv[2]<-lambda*t*alpha/beta;
      mv[3]<-var(Ans3(n,t,lambda,alpha,beta));
      mv[4]<-lambda*t*(alpha+alpha^2)/beta^2;
      return(mv)
      } ##以此输出均值的计算值，理论值，方差的计算值与理论值
    mv(10000,10,5,2,4)
    mv(10000,10,2,1,3)
    mv(10000,10,2,4,0.5)
    mv(50000,10,2,1,3)

## ----eval=TRUE----------------------------------------------------------------
mcBeta <- function(n,x){
  y<-runif(n,0,x);
  mc<-30*x*mean((y*(1-y))^2);
  return(mc)
  } ## n为计算次数,用monte carlo法估计积分

for(i in 1:9){
  cat("x=",i/10,"Estimates=",mcBeta(10000,i/10),"pbeta=",pbeta(i/10,3,3),"\n")
  }  ## 输出该积分从从0.1到0.9的计算值与理论值

## ----eval=TRUE----------------------------------------------------------------
Yray <- function(z,sigma){
  y<-((2*sigma^2)*log(1/(1-z)))^(1/2);
  return(y)
  } ## 生成随机变量的函数

Reduction <- function(n,sigma){
  x<-runif(n);
  x1<-runif(n);
  x2<-runif(n);
  X<-Yray(x,sigma);
  XX<-Yray(1-x,sigma);
  X1<-Yray(x1,sigma);
  X2<-Yray(x2,sigma);
  reduction<-sd((X+XX)/2)/sd((X1+X2)/2);
  return(reduction)
  } ## n为计算的样本个数,sigma为密度函数的参数

## 因此举例减少率为
set.seed(666)
Reduction(10000,1)


## ----eval=TRUE----------------------------------------------------------------
g <- function(x){
  y<-(exp((-1/2)*x^2))*(x^2)/(2*pi)^(1/2);
  return(y)
  } ## 函数g

f1 <- function(x){
  y<-x*exp((1-x^2)/2);
  return(y)
  } ## 函数f1

f2 <- function(x){
  y<-(x-1)*exp((-(x-1)^2)/2);
  return(y)
  } ## 函数f2

h1 <- function(x){
  y<-(1-2*log(1-x))^(1/2);
  return(y)
  } ## 生成pdf为f1的随机变量

h2 <- function(x){
  y<-1+(-2*log(1-x))^(1/2);
  return(y)
  } ## 生成pdf为f2的随机变量

gf1 <- function(n){
  x1<-runif(n);
  X1<-h1(x1);
  y1<-(g(X1))/(f1(X1));
  return(y1)
  } ## 生成n个f1时对应的取值

gf2 <- function(n){
  x2<-runif(n);
  X2<-h2(x2);
  y2<-(g(X2))/(f2(X2));
  return(y2)
  } ## 生成n个f2时对应的取值

## 输出结果
set.seed(1);
xx1 <- gf1(10000);
set.seed(1);
xx2 <- gf2(10000);
cat(" f1:mean=",mean(xx1),"sd=",sd(xx1),"\n","f2:mean=",mean(xx2),"sd=",sd(xx2))
##因此估计次积分的值在0.401,0.402附近，且f1比f2的估计有更小的方差

## ----eval=TRUE----------------------------------------------------------------
set.seed(666)
n<-20;alpha<-0.05
UCL <- replicate(10000, expr ={
  x <- rchisq(n, df=2);
  ucl1<-mean(x)-(sd(x))*qt(1-alpha/2,df=n-1)/(n^0.5);
  ucl2<-mean(x)+(sd(x))*qt(1-alpha/2,df=n-1)/(n^0.5);
  c(ucl1,ucl2)
  }) ## 计算置信区间

cat("覆盖概率为:",mean(UCL[1,]<2&UCL[2,]>2))
## 因此相比6.6中0.773的覆盖率要好

## ----eval=TRUE----------------------------------------------------------------

set.seed(666)
n <- 200;alpha <-0.05;Ierror <- numeric(3);m <- 10000
p<-matrix(0,nrow=m,ncol=3)

for (i in 1:m) {
  xx1 <- rchisq(n,df=1);
  xx2 <- runif(n,min=0,max=2);
  xx3 <- rexp(n,rate=1);
  t1 <- t.test(xx1,alternative = "two.sided", mu = 1);
  t2 <- t.test(xx2,alternative = "two.sided", mu = 1);
  t3 <- t.test(xx3,alternative = "two.sided", mu = 1);
  p[i,1] <- t1$p.value;
  p[i,2] <- t2$p.value;
  p[i,3] <- t3$p.value
  } ##t-test

for(i in 1:3){
  Ierror[i]<-mean(p[,i] < alpha)
  } ## 计算I类错误率

for(i in 1:3){cat("第",i,"个分布的I类错误率为",Ierror[i],"\n")}

## 因此我们可以看到这三个值都与0.05相近,也就说明了t检验对轻微偏离正态分布的检验是稳健的

## ----eval=TRUE----------------------------------------------------------------
set.seed(666)
n <- c(10, 20, 30, 50, 100, 500);m <- 2000;p_reject <- numeric(length(n))
Ans1 <- function(d){
  cv<- qchisq(0.95,d*(d + 1)*(d + 2)/6)
  for(i in 1:length(n)){
    sktests <- numeric(m)
    for(j in 1:m){
      X<-matrix(rnorm(d*n[i]),n[i],d)
      Xs<-scale(X)
      skm<-(tcrossprod(Xs,Xs))^3
      sk<-sum(skm)/(6*n[i])
      sktests[j]<- (abs(sk) >= cv )
    } ## j在m=10000次循环
  p_reject[i] <- mean(sktests)
  } ## i在6个n里循环
  for(ni in 1:length(n)){
    cat("样本数为",n[ni],"时,I类错误概率为",p_reject[ni],"\n")
    } ## 输出
} ## 函数结束

Ans1(2) ##以二维三维的情形为例输出

Ans1(3) 


## ----eval=TRUE----------------------------------------------------------------
alpha <-0.1;n <- 30;m <- 2500
epsilon <- c(seq(0, .14, .01), seq(.15, 1, .05));N <- length(epsilon);p <- numeric(N)
Ans2 <- function(d){
  cv<- qchisq(0.95,d*(d + 1)*(d + 2)/6)
  for (j in 1:N) {
    e <- epsilon[j]
    sktests <- numeric(m)
    for (i in 1:m) {
      s <- sample(c(1, 10), replace = TRUE,size = n*d, prob = c(1-e, e))
      X <- matrix(rnorm(n*d, 0, s),n,d)
      Xs<-scale(X)
      skm<-(tcrossprod(Xs,Xs))^3
      sk<-sum(skm)/(6*n)
      sktests[i] <- (abs(sk) >= cv)
    } ##重复m次
    p[j] <- mean(sktests)
  } ## epsilon结尾
  plot(epsilon, p, type = "b", xlab = bquote(epsilon), ylim = c(0,1))
  abline(h = alpha/2, lty = 3)
} ## 函数结尾

Ans2(2) ## 依次输出二三维的情况

Ans2(3)

## ----eval=TRUE----------------------------------------------------------------

library(boot)
library(MASS)
library(bootstrap)
data(scor, package = "bootstrap")
set.seed(666)
B<-2000
n<-nrow(scor)
theta.b<-numeric(B)
theta.jack<-numeric(n+1)
alpha<-0.05

for(i in 1:(n+1)){
  ev <- eigen(cov(scor[-i,]))
  theta.jack[i]<-(ev$values[1])/sum(ev$values)
  } ## jackknife

theta.hat<-theta.jack[n+1]
theta.jack<-theta.jack[-n-1]
bias.jack<-(n-1)*(theta.hat-mean(theta.jack))
se.jack<-(sd(theta.jack))*(n-1)/(n^0.5)

for(b in 1:B){
  i <- sample(1:n, size = n, replace = TRUE)
  ev <- eigen(cov(scor[i,]))
  theta.b[b]<-(ev$values[1])/sum(ev$values)
  } ## Bootstrap

bias.boot<-mean(theta.b-theta.hat)
se.boot<-sd(theta.b)

round(c(bias_boot=bias.boot,se_boot=se.boot,bias_jack=bias.jack,se_jack=se.jack),5) ##输出bootstrap和jackknife估计下偏差和标准差的估计

perc1 <- as.numeric(quantile(theta.b,alpha/2))
perc2 <- as.numeric(quantile(theta.b,1-alpha/2))

z0 <- qnorm(mean(theta.b<theta.hat))
a <- sum((mean(theta.jack)-theta.jack)^3)/6/sum(((mean(theta.jack)-theta.jack)^2)^1.5)
alpha1 <- pnorm(z0+(z0+qnorm(alpha/2))/(1-a*(z0+qnorm(alpha/2))))
alpha2 <- pnorm(z0+(z0+qnorm(1-alpha/2))/(1-a*(z0+qnorm(1-alpha/2))))

bca1 <- as.numeric(quantile(theta.b,alpha1))
bca2 <- as.numeric(quantile(theta.b,alpha2))

cat("Percentile CI:(",perc1,",",perc2,")") ## 输出percentile置信区间
cat("Bias-corrected and accelerated CI:(",bca1,",",bca2,")") ## 输出BCa置信区间

## ----eval=TRUE----------------------------------------------------------------
library(moments)
set.seed(666)
sk1<-0
sk2<-(1.6)^0.5
n <- 100
B <- 1000
m <- 5000
alpha <- 0.05
ci.norm1<-ci.basic1<-ci.perc1<-matrix(0,m,2)
ci.norm2<-ci.basic2<-ci.perc2<-matrix(0,m,2)
sk1.b<-numeric(B)
sk2.b<-numeric(B)

for(i in 1:m){
  xx1<-rnorm(n)
  xx2<-rchisq(n,5)
  for(s in 1:B){
    sam1 <- sample(1:n, size = n, replace = TRUE)
    sk1.b[s] <- skewness(xx1[sam1])
  } ## bootstrap 1
  for(t in 1:B){
    sam2 <- sample(1:n, size = n, replace = TRUE)
    sk2.b[t] <- skewness(xx2[sam2])
  } ## bootstrap 2
  sk1.hat <- skewness(xx1)
  sk2.hat <- skewness(xx2)
  
  ci.norm1[i,1] <- sk1.hat-qnorm(1-alpha/2)*sd(sk1.b)
  ci.norm1[i,2] <- sk1.hat-qnorm(alpha/2)*sd(sk1.b)
  ci.norm2[i,1] <- sk2.hat-qnorm(1-alpha/2)*sd(sk2.b)
  ci.norm2[i,2] <- sk2.hat-qnorm(alpha/2)*sd(sk2.b)
  
  ci.basic1[i,1] <- 2*sk1.hat-as.numeric(quantile(sk1.b,1-alpha/2))
  ci.basic1[i,2] <- 2*sk1.hat-as.numeric(quantile(sk1.b,alpha/2))
  ci.basic2[i,1] <- 2*sk2.hat-as.numeric(quantile(sk2.b,1-alpha/2))
  ci.basic2[i,2] <- 2*sk2.hat-as.numeric(quantile(sk2.b,alpha/2))
  
  ci.perc1[i,1] <- as.numeric(quantile(sk1.b,alpha/2))
  ci.perc1[i,2] <- as.numeric(quantile(sk1.b,1-alpha/2))
  ci.perc2[i,1] <- as.numeric(quantile(sk2.b,alpha/2))
  ci.perc2[i,2] <- as.numeric(quantile(sk2.b,1-alpha/2))
} ## monte-carlo

norm1.left <- mean(sk1 < ci.norm1[,1] & sk1 < ci.norm1[,2])
norm1.right <- mean(sk1 > ci.norm1[,1] & sk1 > ci.norm1[,2])
norm2.left <- mean(sk2 < ci.norm2[,1] & sk2 < ci.norm2[,2])
norm2.right <- mean(sk2 > ci.norm2[,1] & sk2 > ci.norm2[,2])

basic1.left <- mean(sk1 < ci.basic1[,1] & sk1 < ci.basic1[,2])
basic1.right <- mean(sk1 > ci.basic1[,1] & sk1 > ci.basic1[,2])
basic2.left <- mean(sk2 < ci.basic2[,1] & sk2 < ci.basic2[,2])
basic2.right <- mean(sk2 > ci.basic2[,1] & sk2 > ci.basic2[,2])

perc1.left <- mean(sk1 < ci.perc1[,1] & sk1 < ci.perc1[,2])
perc1.right <- mean(sk1 > ci.perc1[,1] & sk1 > ci.perc1[,2])
perc2.left <- mean(sk2 < ci.perc2[,1] & sk2 < ci.perc2[,2])
perc2.right <- mean(sk2 > ci.perc2[,1] & sk2 > ci.perc2[,2])

c(norm.left=norm1.left,norm.right=norm1.right,basic.left=basic1.left,basic.right=basic1.right,perc.left=perc1.left,perc.right=perc1.right) ## 输出样本服从标准正态分布时的情形

c(norm.left=norm2.left,norm.right=norm2.right,basic.left=basic2.left,basic.right=basic2.right,perc.left=perc2.left,perc.right=perc2.right) ## 输出样本服从自由度为5的卡方分布时的情形

## ----eval=TRUE----------------------------------------------------------------
Ans1 <- function(num,P){
  set.seed(666)
  xxx <- rnorm(num);yyy <- rnorm(num);zzz <- c(xxx,yyy)
  spearman_cor_test = cor.test(xxx, yyy, method = "spearman")
  myl <- rep(0,P);myl0 <- spearman_cor_test$estimate
  for(p in 1:P){
    sam <- sample(1:(2*num), size=num , replace=FALSE)
    xxx1 <- zzz[sam];yyy1 <- zzz[-sam]
    myl[p] <- cor(xxx1, yyy1, method = "spearman")
  } ## permutation test
  p <- mean((abs(c(myl0, myl)) - abs(myl0)) >= 0)
  round(c(p_value_cor=p,p_value_cor.test=spearman_cor_test$p.value),4)
} ## 输入num为两样本的样本元素个数，P为计算次数

Ans1(20,1e3-1) ##以元素样本量为20为例进行输出，通过几组实验认为当计算次数更多后用permutation得到的p值会更接近cor.test得到的p值
suppressWarnings(Ans1(20,1e4-1))
Ans1(20,1e5-1)

## ----eval=TRUE----------------------------------------------------------------

library(RANN)
library(boot)
library(energy)
library(Ball)

set.seed(666)
P<-1e2;alpha<-0.1
nr1 <- nr2<-50;nr <- nr1+nr2;NR <- c(nr1,nr2);sizes<-NR
k<-3;R<-999
p1<-p2<-p3<-rep(0,P)
  
Tn <- function(zz, sam,sizes,k) {
  nr1<-sizes[1];nr2 <- sizes[2];nr<-nr1+nr2
  if(is.vector(zz)) zz <- data.frame(zz,0)
  zz <- zz[sam, ]
  nntests <- nn2(data=zz, k=k+1)
  b1 <- nntests$nn.idx[1:nr1,-1]
  b2 <- nntests$nn.idx[(nr1+1):nr,-1]
  sum1 <- sum(b1 < nr1 + 0.5)
  sum2 <- sum(b2 > nr1 + 0.5)
  tn<-(sum1 +sum2) / (k * nr)
  return(tn)
} ## Tn函数在此结束
  
nn.test <- function(zz,sizes,k){
  objection <- boot(data=zz,statistic=Tn,R=R,sim = "permutation",sizes = sizes,k=k)
  aa <- c(objection$t0,objection$t)
  value_of_p <- mean(aa>=aa[1]) 
  list(statistic=aa[1],p.value=value_of_p)
} ## nn函数在此结束
  
for(j in 1:P){
  xx <- matrix(rnorm(2*nr1,0,1),nr1,2)
  yy <- matrix(rnorm(2*nr2,0,1.38),nr2,2)
  zz <- rbind(xx,yy)
  p1[j] <- nn.test(zz,NR,k)$p.value
  p2[j] <- eqdist.etest(zz,sizes=NR,R=R)$p.value
  p3[j] <- bd.test(xx,yy,num.permutations=999,seed=j*666)$p.value
}
pow1 <- mean(p1<alpha)
pow2 <- mean(p2<alpha)
pow3 <- mean(p3<alpha)
round(c(NN=pow1,Energy=pow2,Ball=pow3),4) ##输出三种方法的power


## ----eval=TRUE----------------------------------------------------------------
set.seed(666)
P<-1e2;alpha<-0.1
nr1 <- nr2<-50;nr <- nr1+nr2;NR <- c(nr1,nr2);sizes<-NR
k<-3;R<-999
p1<-p2<-p3<-rep(0,P)
  
Tn <- function(zz, sam,sizes,k) {
  nr1<-sizes[1];nr2 <- sizes[2];nr<-nr1+nr2
  if(is.vector(zz)) zz <- data.frame(zz,0)
  zz <- zz[sam, ]
  nntests <- nn2(data=zz, k=k+1)
  b1 <- nntests$nn.idx[1:nr1,-1]
  b2 <- nntests$nn.idx[(nr1+1):nr,-1]
  sum1 <- sum(b1 < nr1 + 0.5)
  sum2 <- sum(b2 > nr1 + 0.5)
  tn<-(sum1 +sum2) / (k * nr)
  return(tn)
} ## Tn函数在此结束
  
nn.test <- function(zz,sizes,k){
  objection <- boot(data=zz,statistic=Tn,R=R,sim = "permutation",sizes = sizes,k=k)
  aa <- c(objection$t0,objection$t)
  value_of_p <- mean(aa>=aa[1]) 
  list(statistic=aa[1],p.value=value_of_p)
} ## nn函数在此结束
  
for(j in 1:P){
  xx <- matrix(rnorm(2*nr1,0.25,1),nr1,2)
  yy <- matrix(rnorm(2*nr2,0,1.38),nr2,2)
  zz <- rbind(xx,yy)
  p1[j] <- nn.test(zz,NR,k)$p.value
  p2[j] <- eqdist.etest(zz,sizes=NR,R=R)$p.value
  p3[j] <- bd.test(xx,yy,num.permutations=999,seed=j*666)$p.value
}
pow1 <- mean(p1<alpha)
pow2 <- mean(p2<alpha)
pow3 <- mean(p3<alpha)
round(c(NN=pow1,Energy=pow2,Ball=pow3),4) ##输出三种方法的power

## ----eval=TRUE----------------------------------------------------------------
set.seed(666)
P<-1e2;alpha<-0.1
nr1 <- nr2<-50;nr <- nr1+nr2;NR <- c(nr1,nr2);sizes<-NR
k<-3;R<-999
p1<-p2<-p3<-rep(0,P)
  
Tn <- function(zz, sam,sizes,k) {
  nr1<-sizes[1];nr2 <- sizes[2];nr<-nr1+nr2
  if(is.vector(zz)) zz <- data.frame(zz,0)
  zz <- zz[sam, ]
  nntests <- nn2(data=zz, k=k+1)
  b1 <- nntests$nn.idx[1:nr1,-1]
  b2 <- nntests$nn.idx[(nr1+1):nr,-1]
  sum1 <- sum(b1 < nr1 + 0.5)
  sum2 <- sum(b2 > nr1 + 0.5)
  tn<-(sum1 +sum2) / (k * nr)
  return(tn)
} ## Tn函数在此结束
  
nn.test <- function(zz,sizes,k){
  objection <- boot(data=zz,statistic=Tn,R=R,sim = "permutation",sizes = sizes,k=k)
  aa <- c(objection$t0,objection$t)
  value_of_p <- mean(aa>=aa[1]) 
  list(statistic=aa[1],p.value=value_of_p)
} ## nn函数在此结束
  
for(j in 1:P){
  xx <- matrix(rt(2*nr1,1),nr1,2)
  bi <- as.numeric(runif(2*nr2)>0.5)
  bimodel <- bi*rnorm(2*nr2,-1,2)+(1-bi)*rnorm(2*nr2,1,2)
  yy <- matrix(bimodel,nr2,2)
  zz <- rbind(xx,yy)
  p1[j] <- nn.test(zz,NR,k)$p.value
  p2[j] <- eqdist.etest(zz,sizes=NR,R=R)$p.value
  p3[j] <- bd.test(xx,yy,num.permutations=999,seed=j*666)$p.value
}
pow1 <- mean(p1<alpha)
pow2 <- mean(p2<alpha)
pow3 <- mean(p3<alpha)
round(c(NN=pow1,Energy=pow2,Ball=pow3),4) ##输出三种方法的power

## ----eval=TRUE----------------------------------------------------------------
set.seed(666)
P<-1e2;alpha<-0.1
nr1 <-10; nr2<-90; nr <- nr1+nr2;NR <- c(nr1,nr2);sizes<-NR
k<-3;R<-999
p1<-p2<-p3<-rep(0,P)
  
Tn <- function(zz, sam,sizes,k) {
  nr1<-sizes[1];nr2 <- sizes[2];nr<-nr1+nr2
  if(is.vector(zz)) zz <- data.frame(zz,0)
  zz <- zz[sam, ]
  nntests <- nn2(data=zz, k=k+1)
  b1 <- nntests$nn.idx[1:nr1,-1]
  b2 <- nntests$nn.idx[(nr1+1):nr,-1]
  sum1 <- sum(b1 < nr1 + 0.5)
  sum2 <- sum(b2 > nr1 + 0.5)
  tn<-(sum1 +sum2) / (k * nr)
  return(tn)
} ## Tn函数在此结束
  
nn.test <- function(zz,sizes,k){
  objection <- boot(data=zz,statistic=Tn,R=R,sim = "permutation",sizes = sizes,k=k)
  aa <- c(objection$t0,objection$t)
  value_of_p <- mean(aa>=aa[1]) 
  list(statistic=aa[1],p.value=value_of_p)
} ## nn函数在此结束
  
for(j in 1:P){
  xx <- matrix(rnorm(2*nr1,0,1),nr1,2)
  yy <- matrix(rnorm(2*nr2,0,2),nr2,2)
  zz <- rbind(xx,yy)
  p1[j] <- nn.test(zz,NR,k)$p.value
  p2[j] <- eqdist.etest(zz,sizes=NR,R=R)$p.value
  p3[j] <- bd.test(xx,yy,num.permutations=999,seed=j*666)$p.value
}
pow1 <- mean(p1<alpha)
pow2 <- mean(p2<alpha)
pow3 <- mean(p3<alpha)
round(c(NN=pow1,Energy=pow2,Ball=pow3),4) ##输出三种方法的power

## ----eval=TRUE----------------------------------------------------------------

set.seed(666)
rw.Metropolis <- function(sigma, x0, N){
  xt <- rep(0,N)
  xt[1] <- x0
  k <- 0
  u <- runif(N)
  for (s in 2:N) {
    y <- rnorm(1, xt[s-1], sigma)
        if (u[s] > (dcauchy(y) / dcauchy(xt[s-1]))){
          xt[s] <- xt[s-1]; k <-k+1
        }
        else {
          xt[s] <- y
          } ## 赋值
        } ##游走
  return(list(x=xt, k=k))
} ##函数

N <- 15000
sigma <- 2.5
x0 <- c(-10, -5, 5, 10)
rw1 <- rw.Metropolis(sigma, x0[1], N) 
rw2 <- rw.Metropolis(sigma, x0[2], N)
rw3 <- rw.Metropolis(sigma, x0[3], N)
rw4 <- rw.Metropolis(sigma, x0[4], N)

## print(c(rw1$k, rw2$k, rw3$k, rw4$k)/N) ##输出这四个链的拒绝概率

rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x)
b <- 1001 ##buring sample
a <- seq(0.1,0.9,0.1)
Q <- qcauchy(a)
mc <- rw[b:N,]
Qrw <- apply(mc,2,function(x) quantile(x,a))
QQ <- cbind(Q, Qrw)
colnames(QQ) <- c("Cauchy Quantiles", "chain1", "chain2", "chain3", "chain4") 
print(round(QQ, 3)) ## 理论分位数与样本分位数的比较

erg.mean <- function(x){
  lx <- length(x)
  result <-cumsum(x)/c(1:lx)
  return(result)
  } 

phimc <- apply(mc, 2, erg.mean)
plot(phimc[,1],type="l", xlab='N', ylab="phi",ylim=c(-2,4),main='均值遍历图')
lines(phimc[,2],type="l", col='blue')
lines(phimc[,3],type="l", col='red')
lines(phimc[,4],type="l", col='green')

ww <- function(x){
  result <- sum((x-mean(x))^2)
  return(result)
  }

n <- nrow(phimc)
bn <- n/3*sum((colMeans(phimc)-mean(phimc))^2)
wn <- 1/4/n*sum(apply(phimc, 2, ww))
hat_r <- round(1-1/n+bn/wn/n,3)
cat("Gelman-Rubin统计量为",hat_r,"小于 1.2 ,因此符合监视收敛要求")

## ----eval=TRUE----------------------------------------------------------------

rm(list = ls())
set.seed(666)
alpha <-2; n<- 16; beta <-4

gibbs <- function(x0, y0, N){
  xyg <- matrix(0, N, 2)
  xyg[1,1] <- x0
  xyg[1,2] <- y0
  for(t in 2:N){
    xyg[t,1] <- rbinom(1, n, xyg[t-1,2])
    xyg[t,2] <- rbeta(1, xyg[t,1]+alpha, -xyg[t,1]+n+beta)
  } ## 赋值
  return(xyg)
} ## 函数

N <- 10000
x0 <- c(4, 3, 6, 2)
y0 <- c(0.3, 0.2, 0.5, 0.4)
xy1 <- gibbs(x0[1], y0[1], N)
xy2 <- gibbs(x0[2], y0[2], N)
xy3 <- gibbs(x0[3], y0[3], N)
xy4 <- gibbs(x0[4], y0[4], N)

## f <- function(x){den <- gamma(n+1)/gamma(x[1]+1)/gamma(n-x[1]+1)*(x[2])^(x[1]+alpha-1)*(1-x[2])^(n-x[1]+beta-1);return(den)}

phix <- cbind(xy1[,1], xy2[,1], xy3[,1], xy4[,1])
b <- 1001

phix <- phix[b:N,]

erg.mean <- function(x){
  lx <- length(x)
  result <-cumsum(x)/c(1:lx)
  return(result)
  } 

phimc <- apply(phix, 2, erg.mean)
plot(phimc[,1],type="l", xlab='N', ylab="phi x",main='X 均值遍历图')
lines(phimc[,2],type="l", col='blue')
lines(phimc[,3],type="l", col='red')
lines(phimc[,4],type="l", col='green')

ww <- function(x){
  result <- sum((x-mean(x))^2)
  return(result)
  }

n <- nrow(phimc)
bn1 <- n/3*sum((colMeans(phimc)-mean(phimc))^2)
wn1 <- 1/4/n*sum(apply(phimc, 2, ww))
hat_r1 <- round(1-1/n+bn1/wn1/n,3)
cat("Gelman-Rubin统计量为",hat_r1,"小于 1.2 ,因此X分量符合监视收敛要求")

####

phiy<- cbind(xy1[,2], xy2[,2], xy3[,2], xy4[,2])
b <- 1001

phiy <- phiy[b:N,]

erg.mean <- function(x){
  lx <- length(x)
  result <-cumsum(x)/c(1:lx)
  return(result)
  } 

phimc <- apply(phiy, 2, erg.mean)
plot(phimc[,1],type="l", xlab='N', ylab="phi y",main='Y 均值遍历图')
lines(phimc[,2],type="l", col='blue')
lines(phimc[,3],type="l", col='red')
lines(phimc[,4],type="l", col='green')

ww <- function(x){
  result <- sum((x-mean(x))^2)
  return(result)
  }

n <- nrow(phimc)
bn2 <- n/3*sum((colMeans(phimc)-mean(phimc))^2)
wn2 <- 1/4/n*sum(apply(phimc, 2, ww))
hat_r2 <- round(1-1/n+bn2/wn2/n,3)
cat("Gelman-Rubin统计量为",hat_r2,"小于 1.2 ,因此y分量符合监视收敛要求")

## ----eval=TRUE----------------------------------------------------------------

Term <- function(va, k){
  a <- sum(va^2)
  d <- length(va)
  term <- a*(-a/2)^k/(2*k+1)/(2*k+2)*exp(lgamma((d+1)/2)+lgamma(k+1.5)-lgamma(k+d/2+1)-lgamma(k+1))
  return(term)
} ## 生成第k个的系数

va <- c(1,2) 

Ans1 <- function(k){Term(va,k)} ##以具体的a为例输出系数

cat("此次计算中最后一个能被计算大于0的系数为k=176时，此时系数为:",Ans1(176),"\n而k大于等于177之后系数的值过小而无法计算,k=177时系数为:",Ans1(177))

Ans2 <- function(n){
  x <- matrix(c(0:n),1,n+1)
  s <- sum(apply(x,2,Ans1))
  return(s)
} ## 输出k从0到n时系数之和

N <- 200
sam <- matrix(c(0:N),1,N+1)
x1 <- apply(sam,2,Ans1)
x2 <- apply(sam,2,Ans2)

plot(x1,type="l", xlab='k', ylab="系数",main='系数变化图')
plot(x2,type="l", xlab='N', ylab="和",main='系数之和变化图')
cat("因此总和为:",Ans2(200))

## ----eval=TRUE----------------------------------------------------------------

rm(list = ls())

cc <- function(k,a){((a^2*k)/(k+1-a^2))^0.5}
f <- function(k,a){
  Int <- function(u){
    y <- (1+u^2/k)^((k+1)/(-2))
    return(y)
  } ## 积分内部用的函数
  inc <- integrate(Int,0,cc(k,a))$value
  yy <- 2*inc/(pi*k)^0.5*exp(lgamma((k+1)/2)-lgamma(k/2))
  return(yy)
} ## f(k)=f(k-1)

alpha <- function(k){
  sol <- function(a){
    re <- f(k,a)-f(k-1,a)
    return(re)
  } ##要解的方程
  hh <- 1.4*as.numeric(k==2)+1.7*as.numeric(k==3)+2*as.numeric(k>3) ## 搞个满足uniroot条件的上界
  solution <- uniroot(sol,c(0.1,hh),tol=0.001)$root
  return(solution)
} #求解alpha k,此时解出的解为正解，相应的对于每个k,alpha还有个相反数同样为解,同时0也为解

N <- 100
sam <- matrix(c(2:N),1)
ak <- apply(sam, 2, alpha)
plot(sam,ak,type="l", xlab='k', ylab="alpha k",main='alpha k')


#接下来求解11.4题

alpha2 <- function(k){
  sol2 <- function(a){
    re <- pt(cc(k,a),k)-pt(cc(k-1,a),k-1)
    return(re)
  } 
  hh <- 1.4*as.numeric(k==2)+1.7*as.numeric(k==3)+2*as.numeric(k>3) ## 搞个满足uniroot条件的上界
  solution <- uniroot(sol2,c(0.1,hh),tol=0.001)$root
  return(solution)
}

ak2 <- apply(sam, 2, alpha2)
plot(sam,ak2,type="l", xlab='k', ylab="alpha2 k",main='alpha2 k')

cat("两方法解出来的alpha(k)在程序中计算值差值最大为:",max(abs(ak2-ak)),",非常的小.")

## ----eval=TRUE----------------------------------------------------------------

rm(list = ls())
sam <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
result <- round(sum(sam)/7,4)
cat("lambda收敛于:",result)

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  set.seed(666)
#  trims <- c(0, 0.1, 0.2, 0.5)
#  x <- rcauchy(100)
#  
#  out1 <- vector("list", length(trims))
#  f <- function(trim) mean(x, trim = trim)
#  for (i in seq_along(trims)) {
#  out1[[i]] <- f(trims[i])
#  }

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  set.seed(666)
#  trims <- c(0, 0.1, 0.2, 0.5)
#  x <- rcauchy(100)
#  
#  out2 <- vector("list", length(trims))
#  for (i in seq_along(trims)) {
#  out2[[i]] <- mean(trims[i], x=x)
#  }

## ----echo = TRUE, eval = TRUE-------------------------------------------------

set.seed(666)
attach(mtcars)
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)

MOD <- lapply(formulas, lm)
rsq <- function(mod) summary(mod)$r.squared
r2 <- lapply(MOD, rsq)
cat("四种模型的R方为",unlist(r2))

bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE) 
  mpg[rows] ~ disp[rows]
})

MOD2 <- lapply(bootstraps, lm)
r22 <- lapply(MOD2, rsq)
cat("十次Bootstrap下的的R方为",unlist(r22))

## ----echo = TRUE, eval = TRUE-------------------------------------------------

set.seed(666)

x <- rnorm(1000)
df <- data.frame(replicate(20, sample(x, 100, replace = TRUE)))
sd1 <- vapply(df, sd, numeric(1))

cat("data.frame全为数值的例子下，计算各列的标准差为:",unlist(sd1))


df2 <- data.frame(replicate(20, sample(x, 100, replace = TRUE)),replicate(10, letters[sample(c(1:26), 100, replace = TRUE)]))
resam <- sample(c(1:ncol(df2)))
df2 <- df2[,resam] ## 打乱列的顺序

che <- vapply(df2, is.numeric, logical(1))
df3 <- df2[unlist(che)]
sd2 <- vapply(df3, sd, numeric(1))
cat("data.frame不全为数值的例子下，计算数字各列的标准差为:",unlist(sd2))

## ----echo = TRUE, eval = TRUE-------------------------------------------------

mcsapply <- function(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
  FUN <- match.fun(FUN)
  answer <- parallel::mclapply(X = X, FUN = FUN, ...)
  if (USE.NAMES && is.character(X) && is.null(names(answer)))
      names(answer) <- X
  if (!identical(simplify, FALSE) && length(answer))
      simplify2array(answer, higher = (simplify == "array"))
  else answer
}


## ----echo = TRUE, eval = TRUE-------------------------------------------------

library(Rcpp)
library(microbenchmark)
set.seed(666)

cppFunction('NumericMatrix gibbsC(int N, int thin, int a, int b, int n) {
  NumericMatrix mat(N, 2);
  double x = 0, y = 0.5;
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < thin; j++) {
      x = rbinom(1, n, y)[0];
      y = rbeta(1, x + a, n - x + b)[0];
    }
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  return(mat);
}')

gibbsR <- function(N, thin, a, b, n){
  mat <- matrix(nrow = N, ncol = 2)
  x <- 0
  y <- 0.5
  for (i in 1:N) {
    for (j in 1:thin) {
      x <- rbinom(1, n, y)
      y <- rbeta(1, x + a, n - x + b)
    }
    mat[i, ] <- c(x, y)
  }
  return(mat)
}

a <- 2
b <- 4
n <- 16

gibbR = gibbsR(1000, 10, a, b, n)
gibbC = gibbsC(1000, 10, a, b, n)

xR <- gibbR[,1]
yR <- gibbR[,2]
xC <- gibbC[,1]
yC <- gibbC[,2]
qqplot(xR, yR, main='R Plot')
qqplot(xC, yC, main='C Plot')

ts <- microbenchmark(gibbR = gibbsR(1000, 10, a, b, n),
                     gibbC = gibbsC(1000, 10, a, b, n))

summary(ts)[,c(1,3,5,6)]

