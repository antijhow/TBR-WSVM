library(MASS)
library(WeightSVM)
library(mclust)

beta.est = function(z.h,y.h)
{
  n = length(y.h)
  tol = 0.0001
  if (n==1)
  {
    temp = z.h%*%t(z.h)
    diag(temp) = diag(temp) + tol
    beta.h = solve(temp)%*%z.h%*%y.h
  } else {
    temp = t(z.h)%*%z.h
    diag(temp) = diag(temp) + tol
    beta.h = solve(temp)%*%t(z.h)%*%y.h
  }
  return (beta.h)
}

sq.res = function(z.l,y,set.h,set.g)
{
  n = dim(z.l)[1]
  lm.dim = dim(z.l)[2]+1
  if (length(set.h)<=lm.dim || length(set.g)<=lm.dim)
  {
    theta.est = beta.est(z.l,y)
    y.est = z.l%*%theta.est
    sse = sum((y - y.est[,1])^2)
  } else {
    z.h = z.l[set.h,]
    y.h = y[set.h]
    z.g = z.l[set.g,]
    y.g = y[set.g]
    theta.h.est = beta.est(z.h,y.h)
    theta.g.est = beta.est(z.g,y.g)
    y.est = array(0,n)
    y.est[set.h] = z.h%*%theta.h.est
    y.est[set.g] = z.g%*%theta.g.est
    res = (y-y.est)^2
  }
  return(res)  
}
train.res = function(z.l,y,set.h,set.g) 
{
  n = dim(z.l)[1]
  lm.dim = dim(z.l)[2]+1
  if (length(set.h)<=lm.dim || length(set.g)<=lm.dim)
  {
    theta.est = beta.est(z.l,y)
    y.est = z.l%*%theta.est
    sse = sum((y - y.est[,1])^2)
  } else {
    z.h = z.l[set.h,]
    y.h = y[set.h]
    z.g = z.l[set.g,]
    y.g = y[set.g]
    theta.h.est = beta.est(z.h,y.h)
    theta.g.est = beta.est(z.g,y.g)
    y.est = array(0,n)
    y.est[set.h] = z.h%*%theta.h.est
    y.est[set.g] = z.g%*%theta.g.est
    res = (y-y.est)
  }
  return(res)  
}


total.sse = function(z.l,y,set.h,set.g) #z.b design matrix of threshold funciton with intercept
{
  n = dim(z.l)[1]
  lm.dim = dim(z.l)[2]+1
  if (length(set.h)<=lm.dim || length(set.g)<=lm.dim)
  {
    theta.est = beta.est(z.l,y)
    y.est = z.l%*%theta.est
    sse = sum((y - y.est[,1])^2)
  } else {
    z.h = z.l[set.h,]
    y.h = y[set.h]
    z.g = z.l[set.g,]
    y.g = y[set.g]
    theta.h.est = beta.est(z.h,y.h)
    theta.g.est = beta.est(z.g,y.g)
    y.est = array(0,n)
    y.est[set.h] = z.h%*%theta.h.est
    y.est[set.g] = z.g%*%theta.g.est
    sse = sum((y-y.est)^2)
  }
  return(sse)  
}

pred.mse = function(z.l,y,set.h,set.g,theta.h,theta.g)
{
  n = dim(z.l)[1]
  z.h = z.l[set.h,]
  z.g = z.l[set.g,]
  y.pred[set.h] = z.h%*%theta.h
  y.pred[set.g] = z.g%*%theta.g
  se = (y-y.pred)^2
}


pred.res = function(z.l,y,set.h,set.g,theta.h,theta.g)
{
  n = dim(z.l)[1]
  z.h = z.l[set.h,]
  z.g = z.l[set.g,]
  y.pred[set.h] = z.h%*%theta.h
  y.pred[set.g] = z.g%*%theta.g
  se = y-y.pred
}

M.value = function(z.l,z.b,y,gamma)
{
  n = length(y)  
  state = matrix(1,n,1)
  state[which(0>z.b%*%gamma)]=-1
  set.h = which(state == 1)
  set.g = which(state == -1)
  y.h = y[set.h]
  y.g = y[set.g]
  z.h = z.l[set.h,]
  z.g = z.l[set.g,]
  if(length(y.h)>ncol(z.l)&& length(y.g)>ncol(z.l))
  {
    theta.h.est = beta.est(z.h,y.h)
    theta.g.est = beta.est(z.g,y.g)
    y.h.est = z.h%*%theta.h.est
    y.g.est = z.g%*%theta.g.est
    M = (sum((y.h-y.h.est[,1])^2)+sum((y.g-y.g.est[,1])^2))/n
  } else if(length(y.h)<=ncol(z.l)) {
    theta.g.est = beta.est(z.g,y.g)
    y.g.est = z.g%*%theta.g.est
    M = sum((y.g-y.g.est[,1])^2)/length(y.g)
  } else {
    theta.h.est = beta.est(z.h,y.h)
    y.h.est = z.h%*%theta.h.est
    M = sum((y.h-y.h.est[,1])^2)/length(y.h)
  }
  return(M)
}


MCMC.series = function(z.l,z.b,y,niter,max,min,q,initial) # max and min can be arrays
{
  d = dim(z.b)[2]
  n = length(y)
  gamma = matrix(0,d,niter)
  gamma[q,] = 1
  gamma[,1] = initial/initial[q]
  for (i in 2:niter) 
  {
    currentx = gamma[,i-1]
    cand = mvrnorm(n=1,currentx[-q],diag(d-1))
    temp = (min <= cand) & (cand <= max)
    test = 1
    if ( any(temp==FALSE) ) test = 0
    proposedx = array(0,d)
    proposedx[-q] = cand
    proposedx[q] = 1
    accept  = min(1,test*exp(-M.value(z.l,z.b,y,proposedx))/exp(-M.value(z.l,z.b,y,currentx)))
    if (runif(1) < accept) gamma[,i] = proposedx else gamma[,i] = currentx
  }
  return(gamma)
}

MCMC = function(z.l,z.b,y,niter,niter1,max,min,q,initial)
{
  gamma = MCMC.series(z.l,z.b,y,niter,max,min,q,initial)
  all.mse = apply(gamma,2,M.value,y=y,z.l=z.l,z.b=z.b)
  g.sel = gamma[,which(all.mse==min(all.mse))[1]]
  S = gamma
  d = nrow(S)
  S.sort = apply(S,1,sort)
  delta = 0.1*(S.sort[nrow(S.sort),]- S.sort[1,])
  min = max = array(0,d)
  for (i in 1:d)
  {
    temp = which(S.sort[,i]==g.sel[i])[1]
    if (temp == 1) 
    {
      min[i] = S.sort[temp,i] - delta[i] 
      max[i] = S.sort[temp+1,i] + delta[i]
    } else if (temp == ncol(S)) {
      min[i] = S.sort[temp-1,i] - delta[i]
      max[i] = S.sort[temp,i] + delta[i]
    } else {
      min[i] = S.sort[temp-1,i] - delta[i]
      max[i] = S.sort[temp+1,i] + delta[i]
    }
  }
  min = min[-q]
  max = max[-q]
  MC.gamma = NULL
  for ( k in 1:niter1 )
  {
    gamma = MCMC.series(z.l,z.b,y,niter,max,min,q,g.sel)
    all.mse = apply(gamma,2,M.value,y=y,z.l=z.l,z.b=z.b)
    g.sel = gamma[,which(all.mse==min(all.mse))[1]]
    MC.gamma = rbind(MC.gamma,g.sel)
    
  }	
  gamma.est = apply(MC.gamma,2,mean)
  return (gamma.est)
}




initial.set = function(z.l,y,up=0.7,low=0.3)
{
  n = length(y)
  m.est = mean(y)
  ini.h = which(y-m.est >0)
  ini.g = which(y-m.est <=0)
  z.h = z.l[ini.h,]
  y.h = y[ini.h]
  z.g = z.l[ini.g,]
  y.g = y[ini.g]
  beta.h = beta.est(z.h,y.h)
  beta.g = beta.est(z.g,y.g)	
  count = 0
  test = 0
  while (test ==0)
  {
    beta.h.pocket = beta.h
    beta.g.pocket = beta.g
    sq.h = (y-(z.l%*%beta.h)[1:n,])^2 
    sq.g = (y-(z.l%*%beta.g)[1:n,])^2
    w = sq.h - sq.g 	
    temp.percent = quantile(w, probs = c(low,up))
    ini.h = which(w < temp.percent[1])
    ini.g = which(w >= temp.percent[2])
    z.h = z.l[ini.h,]
    y.h = y[ini.h]
    z.g = z.l[ini.g,]
    y.g = y[ini.g]
    beta.h = beta.est(z.h,y.h)
    beta.g = beta.est(z.g,y.g)
    temp = sqrt(sum((beta.h- beta.h.pocket)^2)+sum((beta.g- beta.g.pocket)^2))
    count = count + 1
    if (temp <=0.0001 || count ==10 ) test = 1
  }
  return (list(initial.h = ini.h,initial.g=ini.g))
  
}

initial.set.kmeans = function(z.l,z.b,y,centers=5)
{
  p = dim(z.b)[2]
  data.group = z.b
  data.group = as.data.frame(data.group)
  kmeans.cluster = kmeans(data.group, centers=centers)
  state = kmeans.cluster[["cluster"]]
  mse.compare = array(Inf,centers)
  mse.compare.max = array(0,centers)
  k.sel = NULL
  for (k in 1:centers)
  {
    if (length(which(state==k))> (p+3)) k.sel = c(k.sel,k)
  }
  for (k in k.sel)
  {
    set.k = which(state == k)
    z.k = z.l[set.k,]
    y.k = y[set.k]
    beta.k = beta.est(z.k,y.k)
    y.k.est = z.k%*%beta.k
    mse.compare[k] = 	mean((y.k - y.k.est[,1])^2)
  }
  k.min = which(mse.compare==min(mse.compare))[1]
  set.h = which(state == k.min)
  k.sel = setdiff(k.sel,k.min)
  z.h = z.l[set.h,]
  y.h = y[set.h]
  beta.h = beta.est(z.h,y.h)
  for (k in k.sel)
  {
    set.k = which(state ==k)
    z.k = z.l[set.k,]
    y.k = y[set.k]
    y.k.est = z.k%*%beta.h
    mse.compare.max[k] = mean((y.k-y.k.est[,1])^2)
  }
  k.max = which(mse.compare.max == max(mse.compare.max))
  set.g = which(state == k.max)
  return (list(initial.h = set.h,initial.g=set.g))
}

threshold.wsvm = function(z.l,z.b,y,initial.h,initial.g,cost=10,upper=50,method="linear")
{
  n = length(y)
  set.h = initial.h
  set.g = initial.g
  loss.min = Inf
  count = 0 
  while (count <= upper)
  {	
    z.h = z.l[set.h,]
    z.g = z.l[set.g,]
    y.h = y[set.h]
    y.g = y[set.g]
    theta.h.est = beta.est(z.h,y.h)
    theta.g.est = beta.est(z.g,y.g)
    y.h.est = z.l%*%theta.h.est
    y.g.est = z.l%*%theta.g.est
    res.sq.h = (y - y.h.est[,1])^2
    res.sq.g = (y - y.g.est[,1])^2
    set.h.svm = which(res.sq.h < res.sq.g)
    state = array(2,n)
    state[set.h.svm] = 1 
    weights = abs(res.sq.h - res.sq.g)
    y.state = factor(state)
    data = data.frame(y = y.state,z.b[,-1])
    bound.wsvm = wsvm(y~.,data=data, weight = weights, kernel = method,cost = cost, scale = FALSE,coef0=1,gamma=1)
    pred.y = predict(bound.wsvm,data)
    pred.y = as.numeric(pred.y)
    set.h = which(pred.y == 1)
    set.g = which(pred.y == 2)	
    loss.temp = total.sse(z.l,y,set.h,set.g)
    test = abs(length(set.h)-length(set.g))
    #if (loss.min > loss.temp)
    #{
    #	loss.min = loss.temp
    #	final.wsvm = bound.wsvm
    #} else test = 1
    if (loss.min > loss.temp && test !=n)
    {
      loss.min = loss.temp
      final.wsvm = bound.wsvm
    } else {
      count = count +1
      temp = initial.set(z.l,y,up=runif(1,0.6,0.9),low=runif(1,0.1,0.4))
      set.h = temp$initial.h
      set.g = temp$initial.g
    }	
    if (loss.min == Inf) count = count -1
  }
  return (list(final.wsvm = final.wsvm))
}



