
#######################################################
# Simulating The bounds of Valuation Distributions in Auctions
# Analysis strategy taken from 
######################################################

setwd("~/GitHubProj/Auctions/ValuationUpperBounds")
#clearing the workspace
#install.packages("expm")
library(pacman)
p_load(MASS, expm) #for matrix sqrt, sqrtm()


#################################
#creating data1
################################

#parameters
#B = 1000   #number of bootstrap repetitions
n = 6     #number of bidders
T.n = 300 #number of auctions per bootstrap repetition
lamda = 0 #probability make a jump bid


#value = rlnorm(n,meanlog = 3, sdlog = 1 )
value.full = matrix(rlnorm(T.n*n,meanlog = 3, sdlog = 1 ),n,T.n)
bid.high.ma = matrix(0,n,T.n)

delta.final = vector(mode = "numeric",length = T.n)
s=1+s
#playing game
for (s in 1:T.n) {
  value = value.full[,s]
  reserve = max(min(value) - 1,1) #must be careful, small draws possible
              #and negative delta creates infinite loops
  standbid = reserve
  bid.high.vec = rep(reserve,n)
  delta = min(1,.05*standbid)
  active = 1:n #set of active players
  out = NULL
  bidder.last = NULL
  count = 0
  
  #single auction game
  while(standbid < (max(value)) ) {
    count = count+1
    if (length(out) < (n-1)) {
      if (length(c(out,bidder.last)) < 1) {
        bidder.next = as.integer(sample(active,1) )
      } else {
        bidder.next = as.integer(sample(active[-c(out,bidder.last)], size = 1))
      }
    } else {
      break
    }
    delta = min(1,.05*standbid)
  
    if (value[bidder.next] > (standbid + delta) ) {
      bidder.last = as.integer(bidder.next)
      temp = rbinom(1,1,lamda)
      if (temp > 0) {
        standbid = runif(1,standbid,value[bidder.next])
        bid.high.vec[bidder.next] = standbid
      } else {
        standbid = standbid + delta
        bid.high.vec[bidder.next] = standbid
      }
    } else {
      if (!(bidder.next %in% out)) {
        out = as.integer(c(out,bidder.next))
      }
    bidder.last = as.integer(bidder.next)
    
    } #  
    if (count > 10000) { 
      browser()
      break
    }
  } #close while loop
  
  #results from one auction
  bid.high.ma[,s] = bid.high.vec
  delta.final[s] = delta
  print(s)
} #close T.n loop

#check data
value.full[,1:10]

write(bid.high.ma, "highbid.txt", ncol = n)



##############################################################
##############################################################
## Generating initial cdf estimates
##############################################################

#first ordering the bids
bid.high.ma = as.matrix(bid.high.ma)
bid.order = apply(bid.high.ma,2,sort)

G.hat.i.n <- function(w,mydat = bid.order) {
  temp = mydat
  temp[bid.order>w] = 0
  temp[bid.order<=w] = 1
  return( (apply(temp,1,sum)/T.n) ) #becomes J by 1 vector
}
  
G.delta <- function(w,mydat = bid.order,delts = delta.final){
  temp = mydat[n,]
  temp[temp>(rep(w,T.n)-delts)] = 0
  return( (sum(temp>0)/T.n) )
}  



phi.func <- function(H,i,n) {
  Left = H*factorial((n-i))*factorial((i-1))/factorial(n)
  if (Left == 0){
    return(0)
  } else {
  out = uniroot( function(phi) Left - as.numeric(integrate(
      function(s) s^(i-1)*(1-s)^(n-i),0,phi)[1]),interval = c(0,1) )
    return( out[1] )
  }
}
  

#creating upper matrix to min over
#function is such that default value of mydat is bid.order
my.upper <- function(w,mydat = bid.order) {
  temp = rep(0,n)
  G.vec = G.hat.i.n(w,mydat)
  for (i in 1:n) {
    temp[i] = phi.func(G.vec[i],i,n)
  } #close for
  return(as.numeric(temp) )
}#close function

#creating lower matrix to max over
#note, if the number of bidding participants is fixed at 6, the
#lower bound is not maxed over and is simply

my.lower <- function(w,mydat = bid.order) {
  return( as.numeric(phi.func(G.delta(w,mydat), (n-1), n) ) )
}

#############################################################
#############################################################
#Raw estimates, without adjustment
#############################################################

#true distribution: lognormal
#plotting true

supp.w = seq( 0, 200, length = T.n )

raw.up = apply(sapply(supp.w,my.upper),2,min)
raw.low = sapply(supp.w,my.lower)
cons.up = apply(sapply(supp.w,my.upper),2,max) 

#plotting
output = plnorm(supp.w,meanlog = 3,sdlog = 1)  #true
plot(supp.w,output, type = 'l',main = paste(
  "Uncorrected HT Estimator: ",T.n," Auctions, ", n, " bidders"))
lines(supp.w,raw.up, col = 'red')
lines(supp.w,raw.low, col = 'blue')
lines(supp.w,cons.up,col = 'green')
legend(x = 'bottomright',c("true","raw upper","raw lower", "conservative upper"),lty = c(1,1),lwd = c(2.5,2.5,2.5), col = c("black","red","blue","green"))


#############################################################
#############################################################
#repeating analysis from Haile and Tamer
#############################################################

mu.func <- function(pt,rho) {
  bottom = sum(exp(pt*rho))
  top = sum(exp(pt*rho)*pt)
  return(top/bottom)
}    
    
#note, this only affects upper bound, since
#this is the only one we can min over

unf.est = sapply(supp.w,my.upper)
rho.vec = c(-500,-100,-10,1)
HT.up = matrix(0,length(rho.vec),T.n)
for (i in 1:length(rho.vec)) {
  HT.up[i,] = apply(unf.est,2,mu.func,rho.vec[i])
}

Rh = length(rho.vec)
col.vec = rainbow((Rh + 3))


plot(supp.w,output, type = 'l',main = paste("Corrected HT Estimator, ",
          T.n, "Auctions, ", n, " bidders, Rho Value in Legend"))
lines(supp.w,raw.up, col = col.vec[1])
lines(supp.w,raw.low, col = col.vec[2])
lines(supp.w,cons.up,col = col.vec[3])
count = 0

for (k in 4:(Rh+3)) {
    count = count + 1 
    lines(supp.w,HT.up[count,],col = col.vec[k])
}

legend(x = 'bottomright',c("true","raw upper","raw lower", "conservative upper", 
       rho.vec),lty = c(1,1),lwd = rep(2.5,16), 
       col = c("black",col.vec[1],col.vec[2],col.vec[3],col.vec[4:(Rh+3)]))


##########################################################################################
##########################################################################################
##########################################################################################
#Using technique from Chernozhukov, Lee, and Rosen
##########################################################################################

#First note that
F.joint <- function(a,b,p) {
  #is order of 1, b is order of 2
  #note prob(x_a<=p, x_b <= p) = prob(max(x_a,x_b) <= p) = prob(x_max{a,b} <= p)
  #so, the following works
  #i.star = max(a,b)
  #G.hat.i.n(p)[i.star]
  
  #but this is faster:
  mean((bid.order[a,] <= p)&(bid.order[b,] <= p))
  
}

 
var.covar <-function(p) {
  out = matrix(0,n,n)
  vec.star = G.hat.i.n(p)
  phi.vec = as.numeric(mapply(phi.func,vec.star,1:n,rep(n,n)))
  multip = matrix(0,n,n)
  
  for(i in 1:n) {
    out[i:n,i] = sapply((i:n),F.joint,i,p) - (vec.star[i])*(vec.star[i:n])   
  }
  
  out = out + t(out)
  diag(out) <- vec.star*(1 - vec.star)
  
  #now for delta method derivative
  diag(multip) <- factorial((n-1) )*factorial((1:n - 1))/
    (factorial(n)*(phi.vec)^(0:(n-1))*(1-phi.vec)^(n-(1:n)) )
  
  out2 = t(multip)%*%out %*% multip
  return(out2)
}

##########################################################################################
##########################################################################################
#note that supp.w must be restricted to areas that return strictly positve estimates
##########################################################################################
#this is a limitation of the CLR methodology as applied here

supp.temp = seq(0,200,length = 1000)
unf.est = sapply(supp.temp,my.upper)
supp.info = which(apply(unf.est, 2, function(s) (sum((0.001 < s)&(s<.9995) )) == n))
supp.info
min = supp.temp[min(supp.info)]
max = supp.temp[max(supp.info)]
supp.w = seq(min,max,length = 200)
out = sapply(supp.w, var.covar)
cov.array = array(data = out, c(n,n,length(supp.w)))


##############################################################
#above was to calculate a covar matrix estimate
#below just plugs in to CLR paper directions
##############################################################

g.v.norm <- function(v,w) {
  return(sqrt( sum(ghat.v(v,w)^2 ) )  ) 
}

ghat.v <- function(v,w) {
  if (abs(det(cov.array[,,w])) < 10^-10) {
    return (rep(0,n))
  }
  return( Re(sqrtm(cov.array[,,w]))[,v] )
}


sn.v <- function(v,w) {
  return( (g.v.norm(v,w)/T.n)      )
}

########################################################################
########################################################################
#step 1
########################################################################
gam.n.ti = 1 - 1/log(T.n) #gamma.n.tilde
R = 1000
Z.r = matrix(rnorm(R*n,0,1),n,R)

#step2 - done in covariance calc above
#step 3 - done, this is ghat.v(v)

############################################################################
############################################################################
#step 4
############################################################################

kn.v <- function(p,w) {
  return( (gam.n.ti - quantile(
        apply(sapply(1:n,ghat.v,w) %*% Z.r,2,max), p)) )
}

p = seq(0,1,length = 100)
length(supp.w)
kn.out = sapply(p,kn.v,supp.w[180])


#from our plot, we see, smaller p has a larfer kn.v
plot(p,kn.out, main = paste("kn.v vs p; ", T.n, " Auctions"))
summary(kn.out)



#calculating theta.hat
#theta.hat = matrix(as.numeric(sapply(supp.w, my.upper) ),n,l.w)



V.hat.func <- function(w,p){
  
  theta.temp = my.upper(w)
  sn = sapply(1:n,sn.v,w)
  which( (theta.temp <= (min(theta.temp + kn.v(p,w)*sn) + 2*kn.v(p,w)*sn) ) )
  
}


############################################################################
############################################################################
#step 5
############################################################################

w = supp.w[180]

kn.v.hat <- function(p,sig,w,set) {
  if (length(set) >1) {
    return(as.numeric(sig - quantile(
      apply( (sapply(1:n,ghat.v,w) %*% Z.r)[set,] ,2,max), p) )     )
  } else {
    return (as.numeric(sig - quantile(
      (sapply(1:n,ghat.v,w) %*% Z.r)[set,], p) )     )
  }
  
}


sig = .95 
p = seq(0,1,length = 100)
set = V.hat.func(supp.w[180],.05)
kn.out.2 = sapply(p,kn.v.hat,sig,supp.w[180],set)

#from our plot, we see, smaller p has a larger kn.v
plot(p,kn.out.2, main = paste("kn.v.hat vs p; ",T.n, " Auctions"))
summary(kn.out.2)

theta.nt <- function(p,sig,w,set) {
  min(my.upper(w) + kn.v.hat(p,sig,w,set)*sapply(1:n,g.v.norm,w)/sqrt(T.n))
  
}
  

  
  
  
############################################################################
############################################################################
#Computing estimator along the support
############################################################################


my.wrapper <- function(w.vec, p1,p2, sig.level) {
  set.list = sapply(w.vec, V.hat.func, p1)
  theta.n0 = mapply(theta.nt,p2,sig.level,w.vec,set.list ) 
}


#may want to resimulate before drawing
#Z.r = matrix(rnorm(R*n,0,1),n,R)
p1 = .04
p2 = .20
sig.level = .95



theta.est = my.wrapper(supp.w,p1,p2,sig.level)


###############################################################]
###############################################################
############################################################
#plotting

output = plnorm(supp.w,meanlog = 3,sdlog = 1)  #true
raw.up = apply(sapply(supp.w,my.upper),2,min)
raw.low = sapply(supp.w,my.lower)
rho.end = -300
unf.est = sapply(supp.w,my.upper)
HT.corrected = apply(unf.est,2,mu.func,rho.end)

plot(supp.w,output, type = 'l',main = paste("CLR vs Haile, Tamer: ",
                                    T.n, "Auctions, ", n, " bidders"))
lines(supp.w,raw.up, col = 'red')
lines(supp.w,raw.low, col = 'blue')
lines(supp.w,theta.est,col = 'green')
lines(supp.w,HT.corrected, col = 'purple')
legend(x = 'bottomright',c("true","raw upper","raw lower", "CLR","HT corrected"),lty = c(1,1),
       lwd = c(2.5,2.5,2.5), col = c("black","red","blue","green","purple"))
                                                                                                             





