set.seed(100)
mat = function(data,T,order)
{
  z.mat = matrix(0,T,order)
  for (i in 1:T) z.mat[i,] = data[i:(order+i-1)]
  z.mat = cbind(1,z.mat)
  return(z.mat)
}
library(dplyr)
library(lubridate)

# data preprocessing
for (i in 2014:2019) {
  file_name = paste0("Total_load_", i, ".csv")
  data = read.csv(file_name)
  data$RowDate = as.Date(data$RowDate, format = '%d/%m/%Y')
  if ("RowTime" %in% colnames(data)) {
    data$RowTime <- NULL
  }
  data = data %>%
    group_by(RowDate) %>%
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))
  data = as.data.frame(data)
  colnames(data) = c("Date", "power")
  assign(paste0("ant_", i), data)
}
ant_2014$Date = as.Date(ant_2014$Date,format='%d/%m/%Y')
ant_2015$Date = as.Date(ant_2015$Date,format='%d/%m/%Y')
ant_2016$Date = as.Date(ant_2016$Date,format='%d/%m/%Y')
ant_2017$Date = as.Date(ant_2017$Date,format='%d/%m/%Y')
ant_2018$Date = as.Date(ant_2018$Date,format='%d/%m/%Y')
ant_2019$Date = as.Date(ant_2019$Date,format='%d/%m/%Y')

data4 = rbind(ant_2014,ant_2015,ant_2016,ant_2017,ant_2018,ant_2019)

mini = min(data4$power)
rang = max(data4$power)-min(data4$power)

data4$power = (data4$power-mini)/rang


p = 14
upper = 15
cost = 50
method="linear"

T = length(data4$power)-p
z.l = mat(data4$power,T,p)
z.b = z.l
y = data4$power[(p+1):(T+p)]


temp = initial.set(z.l,y)
set.h = temp$initial.h
set.g = temp$initial.g
temp = threshold.wsvm(z.l,z.b,y,set.h,set.g,cost = cost,upper=upper,method=method)

pred.y = predict(temp$final.wsvm,z.b[,-1])
pred.y = as.numeric(pred.y)
set.h.final = which(pred.y == 1)
length(set.h.final)
set.g.final = which(pred.y == 2)
length(set.g.final)
mse = total.sse(z.l,y,set.h.final,set.g.final)/T
res = sq.res(z.l,y,set.h.final,set.g.final)
mean(res)

theta.h = beta.est(z.l[set.h.final,],y[set.h.final])
theta.g = beta.est(z.l[set.g.final,],y[set.g.final])	

temp = temp$final.wsvm
coef = temp$coefs[,1]
index.sv = temp$index
bound.sv = z.b[index.sv,]
bound.sv = bound.sv[,-1]
final.est = 0
for (i in 1:length(index.sv))
{
  final.est = final.est + coef[i]*bound.sv[i,]
}
intercept = -temp$rho
final.est = c(intercept,final.est)
theta.b = final.est/sqrt(sum(final.est^2))

theta.h
theta.g
theta.b


n = 100
upper = 15
cost = 50

method="linear"
start = proc.time()
loop = 100

p = 14


acc.rate = array(NA,loop)
balance = array(0,loop)
mse = array(0,loop)
theta.b.matrix = matrix(0,loop,p+1)
for (k in 1:loop)
{
	set.seed(k)
  error = rnorm(n,sd=0.3)
	y.class = array(0,n)
	y = array(0,n+p)
	y[1:p] = rnorm(p,sd=0.3)
	z.b = matrix(0,n,1+p)
	for ( i in (p+1):(n+p) )
	{
		z.b[i-p,] = c(1,y[(i-p):(i-1)])
		if (sum(z.b[i-p,]*theta.b) >= 0)
		{	
			y.class[i-p] = 1
			y[i] = sum(z.b[i-p,]*theta.h) + error[i-p]
		} else {
			y.class[i-p] = 2
			y[i] = sum(z.b[i-p,]*theta.g) + error[i-p]
		}
	}
	true = z.b%*%theta.b
	index.h = which(y.class==1)
	index.g = which(y.class==2)
	balance[k] = length(index.g)/n	
	z.l = z.b
	z.b = z.b[,-1]
	z.b = (z.b - min(z.b))/(max(z.b)-min(z.b))
	z.b = cbind(1,z.b)
	y = y[(p+1):(n+p)]
	#plot(y,main="One Sample Pattern of Experiment VI",xaxt="n",xlab="",xlim=c(1,n),ylab="Response",cex.main=2,cex.lab=1.5)
	#points(y=y[index.g],x=c(index.g),col="blue")
	#points(y=y[index.h],x=c(index.h),col="red",pch=16)
	temp = initial.set(z.l,y)
	set.h = temp$initial.h
	set.g = temp$initial.g
	set.seed(k)
	temp.svm = threshold.wsvm(z.l,z.b,y,set.h,set.g,cost = cost,upper=upper,method=method)
	temp = temp.svm$final.wsvm
	coef = temp$coefs[,1]
	index.sv = temp$index
	bound.sv = z.b[index.sv,]
	bound.sv = bound.sv[,-1]
  
	final.est = 0
	for (i in 1:length(index.sv))
	{
		final.est = final.est + coef[i]*bound.sv[i,]
	}
	intercept = -temp$rho
	final.est = c(intercept,final.est)
	theta.b.matrix[k,] = final.est/sqrt(sum(final.est^2))
	temp = z.b%*%theta.b.matrix[k,]
	test.h = which(temp>=0)
	test.g = which(temp<0)
	acc.rate[k] = length(which(temp*true >=0))/n
	if (acc.rate[k]<0.5) 
	{
	  theta.b.matrix[k,]= -theta.b.matrix[k,]
		acc.rate[k] = 1-acc.rate[k]
	}
	mse[k] = total.sse(z.l,y,test.h,test.g)/n
	print(k)
}


proc.time()-start
apply(theta.b.matrix,2,mean)
apply(theta.b.matrix,2,sd)
theta.b
mean(acc.rate)
sd(acc.rate)
mean(balance)
sd(balance)
mean(mse)
sd(mse)
