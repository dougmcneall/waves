# Functions for hmwave_demo_0.Rmd
# Helper functions, or individual processes within the overall framework.

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

ChooseMaximinNroy = function(n.app, waveobj, nreps){
  # Choose a set of NROY points with the largest minimum
  # distance
  ix.list = vector(mode='list', length = nreps)
  mindist.vec = rep(NA, nreps)
  
  for(i in 1:nreps){
    ix = sample(1:nrow(waveobj$X.nroy), n.app)
    X.cand = waveobj$X.nroy[ix, ]
    ix.list[[i]] = ix
    mindist = min(dist( X.cand))
    mindist.vec[i] = mindist
  }
  ix.chosen = ix.list[[which.max(mindist.vec)]]
  
  return(waveobj$X.nroy[ix.chosen, ])
}

add.nroy.design.points = function(X, Y, Y.target, n.aug,
                                  mins.aug, maxes.aug,
                                  thres = 3, disc.list,
                                  disc.sd.list, obs.sd.list){
  
  # Add NROY design points to a design, using uniform sampling from the
  # entire input space.
  # Inputs
  # X            ...       design matrix (output from maximinLHS)    
  # Y            ...       model output matrix
  # Y.target     ...       Target output, or "observation"
  # n.aug        ...       number of candidate points to augment the lhs
  # mins.aug
  # maxes.aug    ...       Limits on the candidate points
  
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
  X.aug = samp.unif(n.aug, mins = mins.aug, maxes = maxes.aug)
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
  
  # I think this code was wrong
  X.nroy.max = apply(X.nroy, 2, max)
  X.nroy.min = apply(X.nroy, 2, min)
  
  return(list(X.nroy = X.nroy,
              X.aug = X.aug, 
              impl.mat = impl.mat, 
              loo.mse.vec = loo.mse.vec,
              fit.list = fit.list,
              X.nroy.max = X.nroy.max,
              X.nroy.min = X.nroy.min)
  )
}

WithinRange = function(x, maxes, mins){
  # which elements of a vector are between
  # elements of the min and max vectors?
  
  all(x < maxes && x > mins)
}



makeTransparent = function(someColor, alpha=100)
  # Transparent colours for plotting
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

reset <- function() {
  # Allows annotation of graphs, resets axes
  par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
  plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}






