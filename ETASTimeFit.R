
library("PtProcess")


  data<-gsub(" ", "", paste(dataname,".csv"))
  XX<-read.csv(data,head=T)
  XX<-XX[XX$mag>=M0,]
 
  xEtas=c()
  XX$time<-gsub("[Z]","",XX$time)
  XX$time<-gsub("[T]"," ",XX$time)

  

  for(i in 1:nrow(XX)){
    xEtas$time[i]=as.numeric(difftime(XX$time[i],XX$time[1],units="days"))
  }
  xEtas$magnitude<-(XX$mag-M0)##Relative magbitude
  xEtas<-as.data.frame(xEtas)
  xEtas<-xEtas[A:B,]
  xEtas$time<-xEtas$time-xEtas$time[1]
  


#----------------------------------------------------------
#    Define Full Model Object

library("PtProcess")

params <- runif(5,0,1)
params[2]=5*params[2]
initial <- log(params[1:5])



#----------------------------------------------------------
#    Define Full Model Object

##The following is essentially the same as what was in the paper
dmagn_mark <- function(x, data, params){
  #  Gamma distribution
  #  exponential density when params[7]=0
  if (params[7]>0){
    lambda <- etas_gif(data, x[,"time"], params=params[1:5])
    y <- dgamma(x[,"magnitude"], shape=1+sqrt(lambda)*params[7],
                rate=params[6], log=TRUE)
  } else y <- dexp(x[,"magnitude"], rate=params[6], log=TRUE)
  return(y)
}

rmagn_mark <- function(ti, data, params){
  #  Gamma distribution
  #  exponential density when params[7]=0
  if (params[7]>0){
    lambda <- etas_gif(data, ti, params=params[1:5])
    y <- rgamma(1, shape=1+sqrt(lambda)*params[7],
                rate=params[6])
  } else y <- rexp(1, rate=params[6])
  return(list(magnitude=y))
}

TT <- c(0, xEtas$time[length(xEtas$time)]) #Where to integrate over


x <- mpp(data=xEtas, gif=etas_gif,
         mark=NULL,
         params=params, TT=TT,
         gmap=expression(params[1:5]),
         mmap=expression(params[1:5]))


#----------------------------------------------------------
#    Fit Model with Exponential Magnitude Distribution


expmap <- function(y, p){
  #   for exponential distribution
  y$params[1:5] <- exp(p)
  return(y)
}
params <- 4*runif(5,0,1)

initial <- log(params[1:5])
z <- optim(initial, neglogLik, object=x, pmap=expmap,
           control=list(trace=1, maxit=1000))

initial <- z$par
z <- nlm(neglogLik, initial, object=x, pmap=expmap,
         print.level=2, iterlim=500, typsize=initial)

#    write estimates to a new model object x0
x0 <- expmap(x, z$estimate)
Min<-c(z$minimum)
out<-c()
for(i in 1:100){
  params <- 4*runif(5,0,1)
  
  initial <- log(params[1:5])
  z2 <- optim(initial, neglogLik, object=x, pmap=expmap,
              control=list(trace=1, maxit=100))
  
  initial <- z2$par
  z2 <- nlm(neglogLik, initial, object=x, pmap=expmap,
            print.level=2, iterlim=500, typsize=initial)
  Min<-c(Min,z2$minimum)
  
    out[[i]]<-z2
  }

print(logLik(x0))


file=gsub(" ", "", paste("A",A,"B",B,dataname,"TimeETAS.RData"))


save.image(file)


