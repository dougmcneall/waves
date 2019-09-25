ls # hmwave.R
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

# testing Vesrion control 

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

# run the model outside of the function
Y = matrix(NA, nrow = n, ncol = d)
colnames(Y) <- c('y1', 'y2', 'y3')
Y = run.model(X)
Y.target = run.model(X.target)

# The standard deviation of the model output is approximately [5,5,1] and the
# mean is around [12,8,2]

obs.sd.list = list(0.2,0.2,0.1)
disc.list = list(0,0,0)
disc.sd.list = list(0, 0, 0) 
thres = 3

mins.aug = rep(0,4)
maxes.aug = rep(1,4)

# Improvements to the above function from Andrianakis et al. (2015)
# 1) Reduce the range the emulator is fit over, as the waves continue.
# 2) Sample candidate design points from NEAR the existing NROY design points.

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




wave1 = add.nroy.design.points(X = X, 
                               Y = Y, 
                               Y.target = Y.target,
                               n.aug = 50000, 
                               mins.aug = mins.aug,
                               maxes.aug = maxes.aug,
                               thres = 3,
                               disc.list=disc.list,
                               disc.sd.list = disc.sd.list,
                               obs.sd.list = obs.sd.list)

# Exclude any design points outside of the ranges returned as NROY



WithinRange = function(x, maxes, mins){
  # which elements of a vector are between
  # elements of the min and max vectors?
  
  all(x < maxes && x > mins)
}

keep1.ix = which(apply(X, FUN = WithinRange,1,
                      maxes = wave1$X.nroy.max,
                      mins = wave1$X.nroy.min))

X2 = rbind(X[keep1.ix, ] , ChooseMaximinNroy(n.app = n.app, waveobj=wave1, nreps = 10000))
Y2 = run.model(X2)

wave2 = add.nroy.design.points(X = X2,
                               Y = Y2, 
                               Y.target = Y.target,
                               n.aug = 30000,
                               mins.aug = wave1$X.nroy.min,
                               maxes.aug = wave1$X.nroy.max,
                               thres = 3,
                               disc.list = disc.list,
                               disc.sd.list = disc.sd.list,
                               obs.sd.list = obs.sd.list)

keep2.ix = which(apply(X2, FUN = WithinRange,1,
                       maxes = wave2$X.nroy.max,
                       mins = wave2$X.nroy.min))

X3 = rbind(X2[keep2.ix, ], ChooseMaximinNroy(n.app = n.app, waveobj=wave2, nreps = 10000))
Y3 = run.model(X3)

wave3 = add.nroy.design.points(X = X3, 
                               Y = Y3,
                               Y.target = Y.target,
                               n.aug = 20000,
                               mins.aug = wave2$X.nroy.min,
                               maxes.aug = wave2$X.nroy.max,
                               thres = 3,
                               disc.list = disc.list,
                               disc.sd.list = disc.sd.list,
                               obs.sd.list = obs.sd.list)

keep3.ix = which(apply(X3, FUN = WithinRange,1,
                       maxes = wave3$X.nroy.max,
                       mins = wave3$X.nroy.min))

X4 = rbind(X3[keep3.ix, ],
           ChooseMaximinNroy(n.app = n.app, waveobj=wave3, nreps = 10000))

Y4 = run.model(X4)

wave4 = add.nroy.design.points(X = X4, 
                               Y = Y4,
                               Y.target = Y.target,
                               n.aug = 20000,
                               mins.aug = wave3$X.nroy.min,
                               maxes.aug = wave3$X.nroy.max,
                               thres = 3,
                               disc.list = disc.list,
                               disc.sd.list = disc.sd.list,
                               obs.sd.list = obs.sd.list)

keep4.ix = which(apply(X4, FUN = WithinRange,1,
                       maxes = wave4$X.nroy.max,
                       mins = wave4$X.nroy.min))

X5 = rbind(X4[keep4.ix, ],
           ChooseMaximinNroy(n.app = n.app, waveobj=wave4, nreps = 10000))

Y5 = run.model(X5)

wave5 = add.nroy.design.points(X = X5,
                               Y = Y5,
                               Y.target = Y.target,
                               n.aug = 20000, 
                               mins.aug = wave4$X.nroy.min,
                               maxes.aug = wave4$X.nroy.max,
                               thres = 3,
                               disc.list = disc.list,
                               disc.sd.list = disc.sd.list,
                               obs.sd.list = obs.sd.list)

keep5.ix = which(apply(X5, FUN = WithinRange,1,
                       maxes = wave5$X.nroy.max,
                       mins = wave5$X.nroy.min))

X6 = rbind(X5[keep5.ix, ], ChooseMaximinNroy(n.app = n.app, waveobj=wave5, nreps = 10000))
Y6 = run.model(X6)

wave6 = add.nroy.design.points(X = X6,
                               Y = Y6,
                               Y.target = Y.target, 
                               n.aug = 20000, 
                               mins.aug = wave5$X.nroy.min,
                               maxes.aug = wave5$X.nroy.max,
                               thres = 3,
                               disc.list = disc.list,
                               disc.sd.list = disc.sd.list,
                               obs.sd.list = obs.sd.list)


all.nroy = rbind(wave1$X.nroy, wave2$X.nroy, wave3$X.nroy, wave4$X.nroy, wave5$X.nroy, wave6$X.nroy, X.target)

cols = c(rep('black',nrow(wave1$X.nroy)), 
         rep('red',nrow(wave2$X.nroy)),
         rep('purple',nrow(wave3$X.nroy)),
         rep('green',nrow(wave4$X.nroy)),
         rep('grey',nrow(wave5$X.nroy)),
         rep('blue',nrow(wave6$X.nroy)),
          rep('gold', 1)
         
)
cex = c(rep(1, nrow(all.nroy)-1), 2)
  
pairs(all.nroy, xlim = c(0,1), ylim = c(0,1), col = cols, cex = cex)

stop()


# How large is the NROY space? 

waves_list = list(wave1, wave2, wave3, wave4 ,wave5, wave6)

PrintNroyProp = function(obj){
  print(nrow(obj$X.nroy))
}

lapply(waves_list, FUN = PrintNroyProp)


# Minimum Implausibility and Optical Depth plots.


# Could use various strategies:

# Smallest memory (but a large compute time) would be to 
# 1) create a matrix which just has two columns using expand.grid
# 2) swap those into a matrix with the right number of columns,
# with other columns made up 
# OR SIMPLER
# Generate a TAAT matrix
# repeat each row a number of times, and fill the relevant
# columns with random numbers.

TaatRandomDesign = function(X, n, reps, mins, maxes){
  
  # Generate a design that has two-at-a-time features,
  # with (reps) random replications, and minima and
  # maxima for each column decided by mins and maxes
  
  nip <- ncol(X) # number of input parameters
  
  col.ix <- combn(1:nip,2)
  
  out <- NULL
  
 # des.nonrep <- NULL
  
  for(i in 1:ncol(col.ix)){
  
    em.vec1 <- seq(from = mins[col.ix[1,i]], to = maxes[col.ix[1, i]], length.out = n)
    em.vec2 <- seq(from = mins[col.ix[2,i]], to = maxes[col.ix[2, i]], length.out = n)
    
    # This bit generates the two "target" columns (the combinatoral bit)
    des.cols <- as.matrix(expand.grid(em.vec1, em.vec2))
    
    # This repeats the rows of the combinatorial bit a number of times
    des.cols.rep =  des.cols[rep(1:nrow(des.cols), each = reps), ]
    
    # for each row of the des.cols, generate a number of repetitions
    
    # Create a holder with random numbers for this section
    # (should be the same nrows as des.cols.rep)
    holder <- samp.unif(reps*(n^2), mins = mins, maxes = maxes)
    
    # place the des.cols.rep in the correct column
    
    #repholder = matrix(rep(holder,each=reps),nrow=reps)
  
    #holder <- matrix(means, ncol = nip, nrow = nrow(des.cols), byrow = TRUE)
  

    
    mat.part <- holder
    
    colu <- col.ix[,i]
    
    mat.part[, as.matrix(colu[1])] <- des.cols.rep[,1]
    mat.part[, as.matrix(colu[2])] <- des.cols.rep[,2]
    
    #des.nonrep <- rbind(des.nonrep, des.cols)
    
    out <- rbind(out, mat.part)
    
  }
  
  return(list(des = out, ix = col.ix))
  
}
  

n = 21
reps = 20
taat.random = TaatRandomDesign(X6,
                               n = n, 
                               reps = reps,
                               mins = rep(0,4),
                               maxes = rep(1,4)
                               )


pred1 = predict(wave6$fit.list[[1]], newdata = taat.random$des, type = 'UK')

impl1 = impl(em = pred1$mean,
             em.sd = pred1$sd, 
             disc = disc.list[[1]],
             disc.sd = disc.sd.list[[1]],
             obs = Y.target[1],
             obs.sd = obs.sd.list[[1]])

blues = brewer.pal(9, 'Blues')

cplot(x = taat.random$des[1:(reps*n^2), 1 ], y = taat.random$des[1:(reps*n^2), 2 ], z = impl1,
      cols = blues, pch = 20)

# OK, now extract the minimum implausibility at each point.

# A function for doing so:

repmat = matrix(impl1 , nrow = reps) # each column contains the reps.
minvec = apply(repmat,2,min) # now this minvec has the same ordering as col.ix

# The shortcut is just to sample every (reps), rather than outputting another matrix
cplot(x = taat.random$des[seq(from = 1, to = reps*n^2, by = reps), 1 ],
      y = taat.random$des[seq(from = 1, to = reps*n^2, by = reps), 2 ], z = minvec,
           cols = blues, pch = 20, cex = 2)


TaatRandomMulti = function(X, fits, n, reps, mins, maxes,
                           disc.list, disc.sd.list, Y.target, obs.sd.list){
  # Calculates the two-at-a-time implausibility, with other dimensions sampled
  # randomly between mins and maxes
  
  taat.random = TaatRandomDesign(X,
                                 n = n, 
                                 reps = reps,
                                 mins = mins,
                                 maxes = maxes
  )
  
  impmat = matrix(NA, nrow = nrow(taat.random$des), ncol = length(fits))
  predmat.mean = matrix(NA, nrow = nrow(taat.random$des), ncol = length(fits))
  predmat.sd = matrix(NA, nrow = nrow(taat.random$des), ncol = length(fits))
  minmat = matrix(NA, nrow = n*n * ncol(taat.random$ix), ncol = length(fits))
  
  for(i in 1:length(fits)){
    
      pred = predict(fits[[i]], newdata = taat.random$des, type = 'UK')
      
      imp = impl(em = pred$mean,
               em.sd = pred$sd, 
               disc = disc.list[[i]],
               disc.sd = disc.sd.list[[i]],
               obs = Y.target[i],
               obs.sd = obs.sd.list[[i]]
              )
    impmat[, i] = imp
    predmat.mean[, i] = pred$mean
    predmat.sd[, i] = pred$sd
    
    repmat = matrix(imp , nrow = reps) # each column contains the reps.
    minvec = apply(repmat,2,min) # now this minvec has the same ordering as col.ix
    minmat[, i] = minvec
    
  }
  
  # Actually find the min too, and sample
  
  return(list(des = taat.random$des, ix =  taat.random$ix, 
            impmat = impmat, predmat.mean = predmat.mean, predmat.sd = predmat.sd,
            minmat = minmat,
            n = n, reps = reps))
  
}
# Repeat over all the outputs

test = TaatRandomMulti(X = X2, fits = wave2$fit.list,
                       n = 21, 
                       reps = 30,
                       mins = rep(0,4),
                       maxes = rep(1,4),
                       disc.list = disc.list,
                       disc.sd.list = disc.sd.list,
                       Y.target = Y.target,
                       obs.sd.list = obs.sd.list)

reset = function() {
  par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
  plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}

PlotMinImpTaat = function(taat, cols){
  
  # setup the layout matrix
  b = matrix(0, ncol(taat$des) - 1, ncol(taat$des) - 1 )
  b[lower.tri(b)| row(b)==col(b)] <- 1:ncol(taat$ix)
  nf <- layout(b)
  
  npc <- taat$n*taat$n
  
  for(i in 1:ncol(taat$ix)){
    
    # indexes into the original design
    i.ix <- ((i * npc) - (npc - 1)) : (i* npc)
    
    # this indexes into the repeated design
    j.ix = (i.ix*taat$reps)-(taat$reps-1)
    
    x.ix <- taat$ix[1,i]
    y.ix <- taat$ix[2,i]
    
    cplotShort(taat$des[j.ix, x.ix],
               taat$des[j.ix, y.ix],
               z =  taat$minmat[i.ix],
               col = cols,
               pch = 20,
               cex = 2,
               xlab = colnames(X)[x.ix],
               ylab = colnames(X)[y.ix]
               #axes = FALSE
    )
  #  points(X.target[x.ix], X.target[y.ix], col = 'black', bg = 'green', cex = 2, pch = 21)
  }
  
}

PlotMinImpTaat(test, cols = blues)

# Need to set maximum implausibility (or at least rescale all the plots in colour)
# Need to add a colour bar.

reset()
image.plot()



ExtractMinImp = function(taat, n, reps){
  
  # the map (order) for the columns we are looking at is in taat$ix
  # Each block of (reps * n^2) has n reps, so we can sample every reps within that
  
  
  
}





                              # mins = apply(X6,2,min), 
                              # maxes = apply(X6,2,max))

# Build predictions using emulators from the last wave.

# The first reps*n^2 will from col.ix[,1] etc.

# to plot this, we 

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


# New function to add candidate points, efficiently in NROY space.

# 
# 1) generate multivariate normal points around current design points
# 2) Check acceptance rate (prop)
# 3) Repeat until acceptance rate is around 20%




# function output needs to be fraction of NROY points, in order to optimise.
# Inputs need to be the standard deviation (i.e. the diagonal of the covariance matrix),
# the set of inputs and the outputs (or the emulator object) 

GenerateNroyMvn = function(x, n.cand, Sigma, fit.list, Y.target,
                           disc.list,
                           disc.sd.list,
                           obs.sd.list,
                           thres) {
  
  # I can't make the optimisation work - try
  # with a single value of sigma
  
  #Sigma = diag(rep(sigma, 4), 4, 4)
  # takes fit, a DiceKriging fit
  # Generate a candidate set of input points, with the aim of keeping
  # a small number of them.
  # This is pretty fast, as the DiceKriging fit is already done
  X.cand = mvrnorm(n = n.cand, mu = x, Sigma = Sigma)
  # x is the NROY point we have so far
  
  # predicting the output at each design point
  pred.list = lapply(fit.list, FUN = 'predict', newdata = X.cand, type = 'UK')
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
  
  nroy.tf = apply(impl.mat, 1, FUN = all.bt, thres = thres)
  nroy.ix = which(nroy.tf==TRUE)
  X.nroy = X.cand[nroy.ix, ]
  
  prop.nroy = nrow(X.nroy) / n.cand
  
  return(list(prop.nroy = prop.nroy))
  # Find the proportion of NROY
 
}

# first, find the NROY points in the first wave

# Generate mvrnorm
Sigma = diag(0.05, 4, 4)
#sigma = 0.05
x = wave1$X.nroy[1,]

p = GenerateNroyMvn(x=x,
                    n.cand = 1000,
                    Sigma=Sigma,
                    fit.list = wave1$fit.list,
                    Y.target=Y.target,
                    disc.list=disc.list,
                    disc.sd.list=disc.sd.list,
                    obs.sd.list=obs.sd.list,
                    thres=3)

GenerateNroyMvnOptimWrap = function(v, x, n.cand, fit.list, Y.target,
                                    disc.list,
                                    disc.sd.list,
                                    obs.sd.list,
                                    thres){
  
  Sigma = diag(v, length(v), length(v))
  
  p = GenerateNroyMvn(x=x,
                      n.cand = n.cand,
                      Sigma=Sigma,
                      fit.list = fit.list,
                      Y.target=Y.target,
                      disc.list=disc.list,
                      disc.sd.list=disc.sd.list,
                      obs.sd.list=obs.sd.list,
                      thres=3)
  out = p$prop.nroy - 0.2 
  out
}

test = GenerateNroyMvnOptimWrap(v = rep(10,4),
                    x=x,
                    n.cand = 10000,
                    fit.list = wave1$fit.list,
                    #Sigma = Sigma,
                    Y.target=Y.target,
                    disc.list=disc.list,
                    disc.sd.list=disc.sd.list,
                    obs.sd.list=obs.sd.list,
                    thres=3)

# Optimize the variances to get a 20% proportion returned

op = optim(par=rep(10,4), fn = GenerateNroyMvnOptimWrap,     
           x=x,
           n.cand = 10000,
           fit.list = wave1$fit.list,
           Y.target=Y.target,
           disc.list=disc.list,
           disc.sd.list=disc.sd.list,
           obs.sd.list=obs.sd.list,
           thres=3,
           control = list(trace = 100),
           method = "SANN"
           )






AddDesignPoints = function(X, Y, Y.target, n.aug, thres = 3,
                                  disc.list, disc.sd.list, obs.sd.list){
  
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
  
  for(i in 1: length(fit.list)){
    
    
    
    
  }
  
  
  
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





# --------------------------------------------------------------
# Mean effect for sensitivity analysis
# --------------------------------------------------------------


MeanEffectDesign = function(design, n, reps, mins, maxes){

# Design to calculate the mean effect of a parameter, while all other
# parameters are varied across their ranges.

# INPUTS:
# design .... original design (e.g. a latin hypercube or output from expand.grid)
# n ......... number of design points in each dimension
# reps ...... number of times to repeat each observation
# un ........ Should we normalise the matrix to un.mins and un.maxes?
# un.mins, un.maxes  ... min and max scaling of the final design

# OUTPUTS:
# ........... (n x nd) rows, nd columns design matrix, sweeping through parameter space
  
oamat <- NULL

nd <- ncol(design)

# mindes <- rep(0, nd)
# maxdes <- rep(1, nd)

for (j in 1:nd){
  # base matrix of 'best' values to put the sweep values into
  #basemat <- matrix(meandes, nrow = n, ncol = nd , byrow = TRUE)
# need to find some way of applying the mins and maxes
  # basemat <- matrix(runif(n*nd*reps), nrow = n*reps, ncol = nd)
  basemat <- samp.unif(n = n*reps, mins = mins, maxes = maxes)
  
  # use seq and length.out
  vec <- seq(from = mins[j], to = maxes[j], length.out = n)
  repvec <- rep(vec, each = reps)
  basemat[ ,j] <- repvec
  oamat <- rbind(oamat,basemat)  
  }

# Apply normalisation to the matrix AFTER it's generated on the unit cube.
#if(un==TRUE){out <- unnormalize(oamat, un.mins = un.mins, un.maxes = un.maxes)}

#else {out <- oamat}

out <- oamat
out

}

# Is there a difference in sensitivity between wave 1 and wave 6?

wave1.maxes = apply(X, 2, max)
wave1.mins = apply(X, 2, min)

reps = 300
nd = 4
n = 50

# Create the design that will calculate the mean effect
X.me.wave1 = MeanEffectDesign(X, n = n, reps = reps, mins = wave1.mins, maxes = wave1.maxes)

# Predictions list for the mean effect
pred.me.wave1  = lapply(wave1$fit.list, FUN = 'predict', newdata = X.me.wave1, type = 'UK')


MeanEffectCalc = function(X.me,y.me, reps, n){
  # Calculate the Mean effects, using output from MeanEffectsDesign
  # And the emulator-predicted output

  # Inputs
  # X.me ...... design output from MeanEffectDesign()
  # y.me ...... emulated mean output at points in X.me
  # reps ...... Number of repeated points at each subdivision
  #             of the parameters
  # n ......... Number of subdivisions of each parameter

  # Output
  # mean.effect.var is the variance of the mean effect, a
  # summary of the mean impact of input [i] on the output
  # taken across other parameter perturbations
  # mean.effect.mat is a matrix with the mean effect of parameter [i]
  # in each column [i]


  reps.mat = matrix(y.me, ncol = reps, byrow = TRUE)
  
  reps.mean = apply(reps.mat,1, mean)

  mean.effect.mat = matrix(reps.mean, nrow = n, byrow = FALSE)

  mean.effect.var = apply(mean.effect.mat, 2, var)

  
  return(list(mean.effect.var = mean.effect.var,
              mean.effect.mat = mean.effect.mat)
         )
}


# Repeat for each of the outputs
me.wave1.y1 = MeanEffectCalc(X.me = X.me.wave1,
  y.me = pred.me.wave1[[1]]$mean,
  reps = reps,
  n = n
  )

me.wave1.y2 = MeanEffectCalc(X.me = X.me.wave1,
  y.me = pred.me.wave1[[2]]$mean,
  reps = reps,
  n = n
  )

me.wave1.y3 = MeanEffectCalc(X.me = X.me.wave1,
  y.me = pred.me.wave1[[3]]$mean,
  reps = reps,
  n = n
  )

 


dev.new(width = 15, height = 5)
par(mfrow = c(1,3))

ylim = range(me.wave1.y1$mean.effect.mat)

plot(1:10, xlim = c(0,1), ylim = ylim, type = 'n')

for(i in 1:nd){

  inseq = seq(from = wave1.mins[i], wave1.maxes[i], length.out = n)
  
  points(inseq, me.wave1.y1$mean.effect.mat[,i], ylim = ylim,
         type = 'o', col = i)

}


ylim = range(me.wave1.y2$mean.effect.mat)
plot(1:10, xlim = c(0,1), ylim = ylim, type = 'n')

for(i in 1:nd){

  inseq = seq(from = wave1.mins[i], wave1.maxes[i], length.out = n)
  
  points(inseq, me.wave1.y2$mean.effect.mat[,i], ylim = ylim,
         type = 'o', col = i)

}


ylim = range(me.wave1.y3$mean.effect.mat)
plot(1:10, xlim = c(0,1), ylim = ylim, type = 'n')

for(i in 1:nd){

  inseq = seq(from = wave1.mins[i], wave1.maxes[i], length.out = n)
  
  points(inseq, me.wave1.y3$mean.effect.mat[,i], ylim = ylim,
         type = 'o', col = i)

}

# Place all of the "variance of mean effects" together to make a matrix

test = rbind(me.wave1.y1$mean.effect.var,
  me.wave1.y2$mean.effect.var,
  me.wave1.y3$mean.effect.var
  )

# Next do the thing for wave 6, and then get rid of the bits that
# are ruled out.



wave6.maxes = apply(X6, 2, max)
wave6.mins = apply(X6, 2, min)


# Create the design that will calculate the mean effect
X.me.wave6 = MeanEffectDesign(X, n = n, reps = reps, mins = wave6.mins, maxes = wave6.maxes)

# Predictions list for the mean effect
pred.me.wave6  = lapply(wave6$fit.list, FUN = 'predict', newdata = X.me.wave6, type = 'UK')



# Repeat for each of the outputs
me.wave6.y1 = MeanEffectCalc(X.me = X.me.wave6,
  y.me = pred.me.wave6[[1]]$mean,
  reps = reps,
  n = n
  )

me.wave6.y2 = MeanEffectCalc(X.me = X.me.wave6,
  y.me = pred.me.wave6[[2]]$mean,
  reps = reps,
  n = n
  )

me.wave6.y3 = MeanEffectCalc(X.me = X.me.wave6,
  y.me = pred.me.wave6[[3]]$mean,
  reps = reps,
  n = n
  )


dev.new(width = 15, height = 5)
par(mfrow = c(1,3))

ylim = range(me.wave6.y1$mean.effect.mat)

plot(1:10, xlim = c(0,1), ylim = ylim, type = 'n')

for(i in 1:nd){

  inseq = seq(from = wave6.mins[i], wave6.maxes[i], length.out = n)
  
  points(inseq, me.wave6.y1$mean.effect.mat[,i], ylim = ylim,
         type = 'o', col = i)

}


ylim = range(me.wave6.y2$mean.effect.mat)
plot(1:10, xlim = c(0,1), ylim = ylim, type = 'n')

for(i in 1:nd){

  inseq = seq(from = wave6.mins[i], wave6.maxes[i], length.out = n)
  
  points(inseq, me.wave6.y2$mean.effect.mat[,i], ylim = ylim,
         type = 'o', col = i)

}


ylim = range(me.wave6.y3$mean.effect.mat)
plot(1:10, xlim = c(0,1), ylim = ylim, type = 'n')

for(i in 1:nd){

  inseq = seq(from = wave6.mins[i], wave6.maxes[i], length.out = n)
  
  points(inseq, me.wave6.y3$mean.effect.mat[,i], ylim = ylim,
         type = 'o', col = i)

}

# Next: know out anything which is ruled out.

  
stop()






               

  



