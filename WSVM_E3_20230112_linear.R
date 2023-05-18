
n = 100
upper = 1
cost = 50

method="linear"
start = proc.time()
loop = 100
theta.h = c(1,0.7,0.2,-1,-1,1,0)
theta.g = c(0.5,-0.4,0.2,0,1,0,-1)
theta.b = c(-2,0.5,0.3,-1,0.6,0,-1) # intercept, other covariates
theta.true = theta.b/sqrt(sum(theta.b^2))
p = 4
q = 2
theta.b.matrix = matrix(0,loop,p+q+1)

acc.rate = array(NA,loop)
balance = array(0,loop)
mse = array(0,loop)


for (k in 1:loop)
{
	set.seed(k)
	error = rnorm(n)
	temp.ind = matrix(rnorm(n*(p-2)),n,p-2)
	temp.cor = mvrnorm(n,mu=c(0,0),Sigma= matrix(c(1,0.8,0.8,1),2,2))
	temp = cbind(temp.ind,temp.cor)
	y.class = array(0,n)
	y = array(0,n+q)
	y[1:q] = array(0,q)
	z.b = matrix(0,n,1+p+q)
	for ( i in (q+1):(n+q) )
	{
		z.b[i-q,] = c(1,y[(i-q):(i-1)],temp[i-q,])
		if (sum(z.b[i-q,]*theta.b) >= 0)
		{	
			y.class[i-q] = 1
			y[i] = sum(z.b[i-q,]*theta.h) + error[i-q]
		} else {
			y.class[i-q] = 2
			y[i] = sum(z.b[i-q,]*theta.g) + error[i-q]
		}
	}		
	index.h = which(y.class==1)
	index.g = which(y.class==2)
	balance[k] = length(index.g)/n	
	

	z.l = z.b
	y = y[(q+1):(n+q)]
	temp = initial.set(z.l,y)
	set.h = temp$initial.h
	set.g = temp$initial.g
	set.seed(k^2)
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
	mse[k] = total.sse(z.l,z.b,y,test.h,test.g)/n
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


