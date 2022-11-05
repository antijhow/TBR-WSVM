
n = 10000
upper = 20
cost = 10
method= "radial"
start = proc.time()
loop = 1
theta.h = c(1,0.5,1)
theta.g = c(-1,-0.5,1)
p = 2
acc.rates = array(0,loop)
balance = array(0,loop)
pred.acc.rates = array(0,loop)
mse = array(0,loop)
mse.pred = array(0,loop)
for (k in 1:loop)
{
	set.seed(k)
	error = rnorm(n)
	z.b = matrix(runif(n*p,-pi,pi),n,p)
	z.b = cbind(1,z.b)
	y.class = array(2,n)
	y.class[which(sin(z.b[,2]*z.b[,3])>=0)] = 1
	index.h = which(y.class==1)
	index.g = which(y.class==2)
	balance[k] = length(index.g)/n	
	
	y = array(0,n)
	y[index.h] = z.b[index.h,]%*%theta.h + error[index.h]
	y[index.g] = z.b[index.g,]%*%theta.g + error[index.g]
	z.l = z.b
	temp = initial.set.kmeans(z.l,z.b,y,centers=sqrt(n))
	set.h = temp$initial.h
	set.g = temp$initial.g
	set.seed(k^2)
	temp = threshold.wsvm(z.l,z.b,y,set.h,set.g,cost = cost,upper,method=method)
	pred.y = predict(temp$final.wsvm,z.b[,-1])
	pred.y = as.numeric(pred.y)
	set.h.final = which(pred.y == 1)
	set.g.final = which(pred.y == 2)
	mse[k] = total.sse(z.l,z.b,y,set.h.final,set.g.final)/n
	
	temp.acc = (length(intersect(set.h.final, index.h))+length(intersect(set.g.final,index.g)))/n
	acc.rates[k] = max(temp.acc,1-temp.acc)
	z.b.pred =  matrix(runif(n*p,-pi,pi),n,p)
	z.b.pred = cbind(1,z.b.pred)
	y.class.pred = array(2,n)
	y.class.pred[which(sin(z.b.pred[,2]*z.b.pred[,3])>=0)] = 1
	index.h.pred = which(y.class.pred==1)
	index.g.pred = which(y.class.pred==2)
	y.pred = array(0,n)
	error.pred = rnorm(n)
	y.pred[index.h.pred] = z.b.pred[index.h.pred,]%*%theta.h + error.pred[index.h.pred]
	y.pred[index.g.pred] = z.b.pred[index.g.pred,]%*%theta.g + error.pred[index.g.pred]
	z.l.pred = z.b.pred
	pred.y.pred = predict(temp$final.wsvm,z.b.pred[,-1])
	pred.y.pred = as.numeric(pred.y.pred)
	set.h.pred = which(pred.y.pred == 1)
	set.g.pred = which(pred.y.pred == 2)
	mse.pred[k] = total.sse(z.l.pred,z.b.pred,y.pred,set.h.pred,set.g.pred)/n	
	pred.temp.acc = (length(intersect(set.h.pred, index.h.pred))+length(intersect(set.g.pred,index.g.pred)))/n
	pred.acc.rates[k] = max(pred.temp.acc,1-pred.temp.acc)

	print(k)
}
proc.time()-start
mean(balance)
sd(balance)
mean(mse)
sd(mse)
mean(acc.rates)
sd(acc.rates)
mean(mse.pred)
sd(mse.pred)
mean(pred.acc.rates)
sd(pred.acc.rates)
