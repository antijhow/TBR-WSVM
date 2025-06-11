
n = 100
upper = 15
cost = 50
method="linear"
start = proc.time()
loop = 100
theta.h = 20*c(n^{-1/2})
theta.g = c(0)
theta.b = c(-1,1,-1,0) # intercept, s1, s2, s3
theta.true = theta.b/sqrt(sum(theta.b^2))
p = 3
theta.b.matrix = matrix(0,loop,p+1)
acc.rate = array(NA,loop)
balance = array(0,loop)
mse = array(0,loop)
for (k in 1:loop)
{
	set.seed(k)
	error = rnorm(n)
	z.b = matrix(0,n,p)
	z.b[,3] = runif(n)
	z.b[,2] = runif(n)
	z.b[,1] = runif(n,-0.5,0.5)+ 1 + z.b[,2]
	z.b = cbind(1,z.b)
	y.class = array(2,n)
	y.class[which(z.b%*%theta.b>=0)] = 1
	index.h = which(y.class==1)
	index.g = which(y.class==2)
	balance[k] = length(index.g)/n	

	y = array(0,n)
	z.l = matrix(1,n,1)
	y[index.h] = z.l[index.h,]*theta.h + error[index.h]
	y[index.g] = z.l[index.g,]*theta.g + error[index.g]
	#temp = initial.set.kmeans(z.l,z.b,y,centers=sqrt(n))
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
	true = z.b%*%theta.true
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
theta.true
mean(acc.rate)
sd(acc.rate)
mean(balance)
sd(balance)
mean(mse)
sd(mse)
