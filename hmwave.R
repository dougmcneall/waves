# hmwave.R
# Functions for producing history matching waves
# History matching process (for each separate output)

# generate initial design
# generate a target point
# generate a uniform sample (for visulasation and calculating space sizes)

# for each output:
# run the model at the design points
# run the model at the target point
# fit a km to each output
# calculate the implausibility of each point in the uniform sample
# calculate the implausibility of each point in the augmented lhs
# Return the model fit and some fit statistics, both sets of 
# implausibility calculations, and the new design

# Could generate a number of waves and then see how the fit
# statistics change

library(lhs)
library(DiceKriging)
source('https://raw.githubusercontent.com/dougmcneall/packages-git/5f79ffe749f25c6fc39f4f7925e1538d36b7caf1/emtools.R')
source('https://raw.githubusercontent.com/dougmcneall/packages-git/5f79ffe749f25c6fc39f4f7925e1538d36b7caf1/imptools.R')
source('https://raw.githubusercontent.com/dougmcneall/packages-git/5f79ffe749f25c6fc39f4f7925e1538d36b7caf1/vistools.R')

# Test functions from the Virtual Library of Simulation Experiments

# Friedman
friedman = function(x){
  
  stopifnot(length(x)==4)
  
  out = 10*sin(pi*x[1]*x[2]) + 20*(x[3] - 0.5)^2 + 10*x[4] # + 5*x5
  
  out
}

park1 = function(x){
  
  stopifnot(length(x)==4)
  
  out = x[1]/2 * sqrt( 1+ (x[2]+ (x[3]^2) * (x[4] / x[1]^2))) + x[1] + 3*x[4]*exp(1+sin(x[3]))
  
  out
}

park2 = function(x){
  
  stopifnot(length(x)==4)
  
  out = (2/3 * exp(x[1]+x[2])) - (x[4] *sin(x[3])) + x[3]
  
  out
}



# Helper functions, or individual processes within the overall framework.
# run the model

run.model = function(X){
  # Generates multivariate output from a design matrix, using
  # the friedman, park1 and park2 functions.
  # Input is a design matrix with 4 columns on the [0,1] cube
  # Output is a 3 column matrix of model output
  stopifnot(ncol(X)==4)
  n = nrow(X)
  
  Y = matrix(NA, nrow = n, ncol = 3)
  colnames(Y) = c('y1', 'y2', 'y3')
  
  for(i in 1:n){
    
    x = X[i, ]
    
    y1 = friedman(x)
    y2 = park1(x)
    y3 = park2(x)
    
    Y[i, ] = c(y1, y2, y3)
    
  }
  Y
}

create.kmfit.list = function(X, Y){
  # create a list of km fits for multivariate output
  Y.list = as.list(as.data.frame(Y))
  
  fit.list = NULL
  for(i in 1:length(Y.list)){
    
    fit = km(formula = ~.,  design = X, response = Y.list[[i]])
    fit.list = c(fit.list, fit)
  }
  fit.list
}

# Are all the elements of a matrix row below a threshold?
all.bt = function(x, thres) all(x < thres)


# where is below the threshold for all cases?
# initial design
n = 10 # number of initial design points
k = 4  # number of model inputs
d = 3  # number of model outputs

# simplest method is to choose a subset of the discovered NROY
# points to append to the design
n.app = 10 #

X = maximinLHS(n, k = k,dup = 2)
colnames(X) <- c('x1', 'x2', 'x3', 'x4')

X.target = matrix(runif(4), nrow = 1)
n.aug = 50000

# run the model outside of the function
Y = matrix(NA, nrow = n, ncol = d)
colnames(Y) <- c('y1', 'y2', 'y3')
Y = run.model(X)
Y.target = run.model(X.target)

# The standard deviation of the model output is approximately [5,5,1] and the
# mean is around [12,8,2]

obs.sd.list = list(0.5,0.5,0.1)
disc.list = list(0,0,0)
disc.sd.list = list(0, 0, 0) 
thres = 3

add.nroy.design.points = function(X, Y, Y.target, n.aug,thres = 3, disc.list,
                                  disc.sd.list, obs.sd.list){
  
  # Add NROY design points to a design, using uniform sampling from the
  # entire input space.
  # Inputs
  # X            ...       design matrix (output from maximinLHS)    
  # Y            ...       model output matrix
  # Y.target     ...       Target output, or "observation"
  # n.aug        ...       number of candidate points to augment the lhs
  
  # list of fitted km objects, one list element for each output
  fit.list = create.kmfit.list(X=X, Y=Y)
  
  # How good is the fit?
  loo.list = lapply(fit.list, FUN = leaveOneOut.km, type = 'UK', trend.reestim = TRUE)
  
  loo.mse.vec = rep(NA, length(loo.list))
  for(i in 1:length(loo.list)){
    
    loo.mse = mean((loo.list[[i]]$mean - Y.target[,i])^2)
    loo.mse.vec[i] = loo.mse
  }
  
  # create a new set of candidate points
  #X.aug = augmentLHS(X, n.aug) # this will add to the current design.
  X.aug = samp.unif(n.aug, mins = rep(0,4), maxes = rep(1,4))
  colnames(X.aug) = colnames(X)
  
  # predicting the output at each design point
  pred.list = lapply(fit.list, FUN = 'predict', newdata = X.aug, type = 'UK')
  # calculate the implausibility at each design point
  
  impl.mat = NULL
  for(i in 1:length(pred.list)){
    
    pred.impl = impl(em = pred.list[[i]]$mean,
                     em.sd = pred.list[[i]]$sd,
                     disc = disc.list[[i]],
                     disc.sd = disc.sd.list[[i]],
                     obs = Y.target[[i]],
                     obs.sd = obs.sd.list[[i]])
    
    impl.mat = cbind(impl.mat, pred.impl)
  }
  
  # Which of the candidte design points are NROY?
  
  # create a matrix of the implausibility measures
  # find the indices of the matrix where all are below the threshold.
  nroy.tf = apply(impl.mat, 1, FUN = all.bt, thres = thres)
  nroy.ix = which(nroy.tf==TRUE)
  X.nroy = X.aug[nroy.ix, ]
  
  return(list(X.nroy = X.nroy,
              X.aug = X.aug, 
              impl.mat = impl.mat, 
              loo.mse.vec = loo.mse.vec,
              fit.list = fit.list
  )
  )
}




test1 = add.nroy.design.points(X = X, Y = Y, Y.target = Y.target, n.aug = n.aug, thres = 3,
                               disc.list=disc.list,
                               disc.sd.list = disc.sd.list,
                               obs.sd.list = obs.sd.list)

X2 = rbind(X,test1$X.nroy[1:n.app, ])
Y2 = run.model(X2)

test2 = add.nroy.design.points(X = X2, Y = Y2, Y.target = Y.target, n.aug=n.aug, thres = 3,
                               disc.list = disc.list,
                               disc.sd.list = disc.sd.list,
                               obs.sd.list = obs.sd.list)

X3 = rbind(X2, test2$X.nroy[1:n.app, ])
Y3 = run.model(X3)

test3 = add.nroy.design.points(X = X3, Y = Y3, Y.target = Y.target, n.aug=n.aug, thres = 3,
                               disc.list = disc.list,
                               disc.sd.list = disc.sd.list,
                               obs.sd.list = obs.sd.list)

X4 = rbind(X3, test3$X.nroy[1:n.app, ])
Y4 = run.model(X4)

test4 = add.nroy.design.points(X = X4, Y = Y4, Y.target = Y.target, n.aug=n.aug, thres = 3,
                               disc.list = disc.list,
                               disc.sd.list = disc.sd.list,
                               obs.sd.list = obs.sd.list)

X5 = rbind(X4, test4$X.nroy[1:n.app, ])
Y5 = run.model(X5)

test5 = add.nroy.design.points(X = X5, Y = Y5, Y.target = Y.target, n.aug=n.aug, thres = 3,
                               disc.list = disc.list,
                               disc.sd.list = disc.sd.list,
                               obs.sd.list = obs.sd.list)

X6 = rbind(X5, test5$X.nroy[1:n.app, ])
Y6 = run.model(X6)

test6 = add.nroy.design.points(X = X6, Y = Y6, Y.target = Y.target, n.aug=n.aug, thres = 3,
                               disc.list = disc.list,
                               disc.sd.list = disc.sd.list,
                               obs.sd.list = obs.sd.list)


all.nroy = rbind(test1$X.nroy, test2$X.nroy, test3$X.nroy, test4$X.nroy, test5$X.nroy, test6$X.nroy, X.target)

cols = c(rep('black',nrow(test1$X.nroy)), 
         rep('red',nrow(test2$X.nroy)),
         rep('purple',nrow(test3$X.nroy)),
         rep('green',nrow(test4$X.nroy)),
         rep('grey',nrow(test5$X.nroy)),
         rep('blue',nrow(test6$X.nroy)),
          rep('gold', 1)
         
)
cex = c(rep(1, nrow(all.nroy)-1), 2)
  
pairs(all.nroy, xlim = c(0,1), ylim = c(0,1), col = cols, cex = cex)

# How large is the NROY space? 

waves_list = list(test1, test2, test3, test4 ,test5, test6)

PrintNroyProp = function(obj){
  print(nrow(obj$X.nroy))
}

lapply(waves_list, FUN = PrintNroyProp)



# Improvements to the above function from Andrianakis et al. (2015)
# 1) Reduce the range the emulator is fit over, as the waves continue.
# 2) Sample candidate design points from NEAR the existing NROY design points.

# from Andrianakis et al. (2015): 
# Suppose that in wave we have a number of non-implausible points. 
# For each of these, we draw samples from a variate normal distribution 
# that is centered on the value of the generating point. The wave 
#implausibility is then evaluated on the new samples and the variance 
# of the normal distribution is selected so that a small percentage of them
# (around 20%) are non-implausible. The low acceptance rates should ensure 
# that the new samples are sufficiently different from the old ones. This 
# method can efficiently generate an adequate number of data points that 
# can be used in subsequent waves.
# A subset of the non-implausible samples drawn are then used to run the 
# simulator and repeat another wave of history matching.





