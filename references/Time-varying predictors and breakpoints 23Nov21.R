# This investigates approximations to the LRT

library(mvtnorm)
library(lattice)
library(raster)
library(Kendall)

source("remote_sensing_tools_1Feb22.R")

###############################
# single time series
###############################
par(mfrow=c(2,2))

AR.model <- T

n.start <- 10
n.obs <- 30

b0 <- 0
be <- 0
bd <- .2
f <- .1
c <- 1
se <- 1
sd <- 1

df <- data.frame(t = 1:n.obs, x=NA, x.lag=NA, y=NA)

e <- 0
d <- 0
y <- 0
x <- 0
for(t in 1:(n.start+n.obs)){
	ee <- rnorm(n=1, sd=se)
	ed <- rnorm(n=1, sd=sd)

	if(t > n.start) df$x.lag[t-n.start] <- x

	e <- be * e + ee
	x.new <- f*t + e

	if(AR.model){
		d <- bd * d + ed
		y <- b0 + c*x.new + d
	}else{
		y <- b0 + bd * y + c*x + ed
	}

	x <- x.new

	if(t > n.start){
		df$x[t-n.start] <- x
		df$y[t-n.start] <- y
	}
}
if(AR.model){
	main=paste0("AR model with bd = ",bd)
}else{
	main=paste0("CLS model with bd = ",bd)
}
plot(y ~ t, data=df, typ="l", ylim=c(min(df[,2:4]), max(df[,2:4])), main=main)
lines(x ~ t, data=df, col="blue")

AR_REML_fit(X = matrix(df$y, nrow=1), t.scale = df$t)$c
AR_REML_fit(X = matrix(df$x, nrow=1), t.scale = df$t)$c

AR_REML_fit(X = matrix(df$y, nrow=1), t.scale = df$x)$c
CLS_fit(X = matrix(df$y, nrow=1), t.scale = df$x)$c
CLS_fit(X = matrix(df$y, nrow=1), t.scale = df$x.lag)$c
LS_fit(X = matrix(df$y, nrow=1), t.scale = df$x)$c

#CLS_fit_U(X = matrix(df$y, nrow=1), U = matrix(df$x, nrow=1))



###############################
# repeats
###############################
AR.model <- T

n.start <- 10
n.obs <- 30

be <- 0
bd <- .8

b0 <- 0
f <- .1
c <- 1
se <- 1
sd <- 1

nreps <- 200
w <- data.frame(1:nreps)
for(i.rep in 1:nreps){

	df <- data.frame(t = 1:n.obs, x=NA, x.lag=NA, y=NA)

	e <- 0
	d <- 0
	y <- 0
	x <- 0
	for(t in 1:(n.start+n.obs)){
		ee <- rnorm(n=1, sd=se)
		ed <- rnorm(n=1, sd=sd)

		if(t > n.start) df$x.lag[t-n.start] <- x

		e <- be * e + ee
		x.new <- f*t + e

		if(AR.model){
			d <- bd * d + ed
			y <- b0 + c*x.new + d
		}else{
			y <- b0 + bd * y + c*x + ed
		}

		x <- x.new

		if(t > n.start){
			df$x[t-n.start] <- x
			df$y[t-n.start] <- y
		}
	}
	w$c.y.x[i.rep] <- AR_REML_fit(X = matrix(df$y, nrow=1), t.scale = df$x)$c
	w$c.x.t[i.rep] <- AR_REML_fit(X = matrix(df$x, nrow=1), t.scale = df$t)$c
	w$c.y.t[i.rep] <- AR_REML_fit(X = matrix(df$y, nrow=1), t.scale = df$t)$c

	w$c.y.x.cls[i.rep] <- CLS_fit(X = matrix(df$y, nrow=1), t.scale = df$x)$c
	w$c.y.x.cls.lag[i.rep] <- CLS_fit(X = matrix(df$y, nrow=1), t.scale = df$x.lag)$c

	w$c.y.x.ls[i.rep] <- LS_fit(X = matrix(df$y, nrow=1), t.scale = df$x)$c
}

summary(lm(c.y.t ~ c.x.t, data=w))$coef
mean(w$c.y.x)
mean(w$c.y.x.cls)
mean(w$c.y.x.cls.lag)
mean(w$c.y.x.ls)

par(mfrow=c(2,2))
plot(c.y.t ~ c.x.t, data=w)
hist(w$c.y.x)
#hist(w$c.y.x.cls)
hist(w$c.y.x.cls.lag)
hist(w$c.y.x.ls)


###############################
# function
###############################
independent <- function(AR.model, bd, nreps){

	n.start <- 10
	n.obs <- 30

	be <- 0

	b0 <- 0
	f <- .1
	c <- 1
	se <- 1
	sd <- 1

	w <- data.frame(1:nreps)
	for(i.rep in 1:nreps){

		df <- data.frame(t = 1:n.obs, x=NA, x.lag=NA, y=NA)

		e <- 0
		d <- 0
		y <- 0
		x <- 0
		for(t in 1:(n.start+n.obs)){
			ee <- rnorm(n=1, sd=se)
			ed <- rnorm(n=1, sd=sd)

			if(t > n.start) df$x.lag[t-n.start] <- x

			e <- be * e + ee
			x.new <- f*t + e

			if(AR.model){
				d <- bd * d + ed
				y <- b0 + c*x.new + d
			}else{
				y <- b0 + bd * y + c*x + ed
			}

			x <- x.new

			if(t > n.start){
				df$x[t-n.start] <- x
				df$y[t-n.start] <- y
			}
		}

		w$c.x.t[i.rep] <- AR_REML_fit(X = matrix(df$x, nrow=1), t.scale = df$t)$c
		w$c.y.t[i.rep] <- AR_REML_fit(X = matrix(df$y, nrow=1), t.scale = df$t)$c

		w$c.y.x[i.rep] <- AR_REML_fit(X = matrix(df$y, nrow=1), t.scale = df$x)$c

		w$c.y.x.cls[i.rep] <- CLS_fit(X = matrix(df$y, nrow=1), t.scale = df$x)$c
		w$c.y.x.cls.lag[i.rep] <- CLS_fit(X = matrix(df$y, nrow=1), t.scale = df$x.lag)$c

		w$c.y.x.ls[i.rep] <- LS_fit(X = matrix(df$y, nrow=1), t.scale = df$x)$c

		w$c.y.x.MSE[i.rep] <- AR_REML_fit(X = matrix(df$y, nrow=1), t.scale = df$x)$MSE
		w$c.y.x.cls.MSE[i.rep] <- CLS_fit(X = matrix(df$y, nrow=1), t.scale = df$x)$MSE
		w$c.y.x.cls.lag.MSE[i.rep] <- CLS_fit(X = matrix(df$y, nrow=1), t.scale = df$x.lag)$MSE
		w$c.y.x.ls.MSE[i.rep] <- LS_fit(X = matrix(df$y, nrow=1), t.scale = df$x)$MSE
	}

	return(c(slope = summary(lm(c.y.t ~ c.x.t, data=w))$coef[2,1], c.AR = mean(w$c.y.x), c.CLS = mean(w$c.y.x.cls), c.CLS.lag = mean(w$c.y.x.cls.lag), c.LS = mean(w$c.y.x.ls), c.AR.MSE = mean(w$c.y.x.MSE), c.CLS.MSE = mean(w$c.y.x.cls.MSE), c.CLS.lag.MSE = mean(w$c.y.x.cls.lag.MSE), c.LS.MSE = mean(w$c.y.x.ls.MSE)))
}

output <- NULL
for(AR.model in c(T,F)) for(bd in c(0,.8)) output <- rbind(output, c(AR.model = AR.model, bd = bd, independent(AR.model, bd, 1000)))
output

     # AR.model  bd    slope       c.AR     c.CLS  c.CLS.lag      c.LS  c.AR.MSE c.CLS.MSE c.CLS.lag.MSE  c.LS.MSE
# [1,]        1 0.0 1.031878  0.9988291 1.0017718  0.4428479 0.9972054 0.9737077  1.009838      2.450332 0.9949481
# [2,]        1 0.8 1.111057  1.0084976 0.8662229 -0.3404192 1.0191618 0.9721931  1.349990      2.444274 1.3537002
# [3,]        0 0.0 0.875013  0.3302709 0.3168172  0.9932594 0.3891311 2.2904810  2.427711      1.001915 1.5557851
# [4,]        0 0.8 1.929508 -0.4598284 0.1399467  1.0189094 2.0758840 2.0502969  2.259988      1.004117 3.8630313


write.table(output, file="output.csv", sep=",", row.names=F)

###############################
# breakpoints
###############################
library(bfast)

AR.model <- T

n.start <- 10
n.obs <- 30

be <- .8
bd <- .8

b0 <- 0.5
f <- .1
c <- 1
se <- 1
sd <- 1

df <- data.frame(t = 1:n.obs, x=NA, x.lag=NA, y=NA)

e <- 0
d <- 0
y <- 0
x <- 0
for(t in 1:(n.start+n.obs)){
	ee <- rnorm(n=1, sd=se)
	ed <- rnorm(n=1, sd=sd)

	if(t > n.start) df$x.lag[t-n.start] <- x

	e <- be * e + ee
	x.new <- f*t + e

	if(AR.model){
		d <- bd * d + ed
		y <- b0 + c*x.new + d
	}else{
		y <- b0 + bd * y + c*x + ed
	}

	x <- x.new

	if(t > n.start){
		df$x[t-n.start] <- x
		df$y[t-n.start] <- y
	}
}
out <- ts(array(df$y),array(df$t), frequency=2)
mod01 <- bfast01(out, formula = response ~ trend, lag = NULL)
mod <- bfast(out)
mod

#plot(mod)
if(AR.model){
	main=paste0("AR model with bx = ",be, " and by = ",bd)
}else{
	main=paste0("CLS model with bx = ",be, " and by = ",bd)
}
par(mfrow=c(2,1), mai=c(1,1,.5,.1))
plot(y ~ t, data=df, typ="l", ylim=c(min(df[,2:4]), max(df[,2:4])), main=main, xlab="Time")
lines(x ~ t, data=df, col="blue")
plot(mod01)


plot(mod)


##################
# repeats
AR.model <- F

n.start <- 10
n.obs <- 30

be <- .8
bd <- .8

b0 <- 0
f <- .1
c <- 1
se <- 1
sd <- 1

nreps <- 1000
w <- data.frame(repp = 1:nreps)
for(i.rep in 1:nreps){

	df <- data.frame(t = 1:n.obs, x=NA, x.lag=NA, y=NA)

	e <- 0
	d <- 0
	y <- 0
	x <- 0
	for(t in 1:(n.start+n.obs)){
		ee <- rnorm(n=1, sd=se)
		ed <- rnorm(n=1, sd=sd)

		if(t > n.start) df$x.lag[t-n.start] <- x

		e <- be * e + ee
		x.new <- f*t + e

		if(AR.model){
			d <- bd * d + ed
			y <- b0 + c*x.new + d
		}else{
			y <- b0 + bd * y + c*x + ed
		}

		x <- x.new

		if(t > n.start){
			df$x[t-n.start] <- x
			df$y[t-n.start] <- y
		}
	}
	out <- ts(array(df$y),array(df$t), frequency=10)
	mod <- bfast01(out, formula = response ~ trend, lag = 1)
	w$reject.lag1[i.rep] <- mod$test

	mod <- bfast01(out, formula = response ~ trend, lag = NULL)
	w$reject.lagNULL[i.rep] <- mod$test

	mod <- bfast(out)
	if(mod$output[[1]]$Vt.bp == 0){
		w$reject.bfast[i.rep] <- 0
	}else{
		w$reject.bfast[i.rep] <- 1
	}
}
i.rep
w <- w[1:i.rep,]
colMeans(w)

# AR with be=.8, bd=.8
         # repp    reject.lag1 reject.lagNULL   reject.bfast
       # 500.500          0.148          0.141          0.097

# CLS with be=.8, bd=.8
          # repp    reject.lag1 reject.lagNULL   reject.bfast
       # 500.500          0.262          0.261          0.212

