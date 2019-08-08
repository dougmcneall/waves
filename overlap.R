# overlap.R
# Exploring using measures of overlap of NROY input space to find model discrepancy.


# Create a simple model of three outputs, with four input variables. In the
# simplest case, all outputs are determined entirely by the four input variables.
# There aren't any feedbacks between the outputs.

# The idea is to add a model discrepancy to one of the outputs, and then see if we can
# identify which output has the discrepancy using history matching, and a measure
# of overlap of the input spaces.
# 
# We should be able to spot the "odd input space out", if the discrepancy is large enough.
# The discrepancy will cause a different input space to be suggested as NROY than the two
# outputs with no discrepancy.

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
  
x = runif(4, min =0, max = 1)

friedman(x)
park1(x)
park2(x)


# Find a latin hypercube sample of x and build a response function (emulator)
library(lhs)

n = 300 # number of design points
k = 4  # number of model inputs
d = 3  # number of model outputs

X = maximinLHS(n, k = k,dup = 2)
colnames(X) <- c('x1', 'x2', 'x3', 'x4')

Y = matrix(NA, nrow = n, ncol = d)
colnames(Y) <- c('y1', 'y2', 'y3')

# Run the "model" at the points in the design
for(i in 1:n){
  
  x = X[i, ]
  
  y1 = friedman(x)
  y2 = park1(x)
  y3 = park2(x)
  
  Y[i, ] = c(y1, y2, y3)
  
}

# Build an emulator of the response surface
library(DiceKriging)

fit.y1 = km(~., design = X, response = Y[, 1])
fit.y2 = km(~., design = X, response = Y[, 2])
fit.y3 = km(~., design = X, response = Y[, 3])


# Sample a point in input space, but don't look where it is
# (we'll try and reconstruct it from the oputput)

x.test = runif(4)
y1.test = friedman(x.test)
y2.test = park1(x.test)
y3.test = park2(x.test)
# load some history matching and emulation tools

source('https://raw.githubusercontent.com/dougmcneall/packages-git/5f79ffe749f25c6fc39f4f7925e1538d36b7caf1/emtools.R')
source('https://raw.githubusercontent.com/dougmcneall/packages-git/5f79ffe749f25c6fc39f4f7925e1538d36b7caf1/imptools.R')
source('https://raw.githubusercontent.com/dougmcneall/packages-git/5f79ffe749f25c6fc39f4f7925e1538d36b7caf1/vistools.R')


# Find the implausibility of a large sample from the emulator.

# sample uniformly from input space
n.unif = 1000
X.unif = samp.unif(n.unif, mins = rep(0,k), maxes = rep(1,k))
colnames(X.unif) <- c('x1', 'x2', 'x3', 'x4')

# predict the model output at each sample point
y1.pred = predict(fit.y1, newdata = X.unif, type = 'UK')
y2.pred = predict(fit.y2, newdata = X.unif, type = 'UK')
y3.pred = predict(fit.y3, newdata = X.unif, type = 'UK')

disc = 0
obs.sd = 0
disc.sd = 1
thres = 3


y1.impl = impl(em = y1.pred$mean, em.sd = y1.pred$sd,
                  disc = disc, obs = y1.test, disc.sd = disc.sd, obs.sd = obs.sd)

y2.impl = impl(em = y2.pred$mean, em.sd = y2.pred$sd,
                disc = disc, obs = y2.test, disc.sd = disc.sd, obs.sd = obs.sd)

y3.impl = impl(em = y3.pred$mean, em.sd = y3.pred$sd,
                disc = disc, obs = y3.test, disc.sd = disc.sd, obs.sd = obs.sd)

impl.mat = cbind(y1.impl, y2.impl, y3.impl)

max.impl = apply(impl.mat, 1, max)

# Find places where everything is below 3

bt.ix = which(max.impl<thres)

X.nroy = X.unif[bt.ix, ]


pairs(X.nroy)

# Does the target input meet the NROY criteria?
y1.test.pred = predict(fit.y1, newdata = matrix(x.test, nrow = 1), type = 'UK')
y2.test.pred = predict(fit.y2, newdata = matrix(x.test, nrow = 1), type = 'UK')
y3.test.pred = predict(fit.y3, newdata = matrix(x.test, nrow = 1), type = 'UK')

y1.test.impl = impl(em = y1.test.pred$mean, em.sd = y1.test.pred$sd,
               disc = disc, obs = y1.test, disc.sd = disc.sd, obs.sd = obs.sd)

y2.test.impl = impl(em = y2.test.pred$mean, em.sd = y2.test.pred$sd,
                    disc = disc, obs = y2.test, disc.sd = disc.sd, obs.sd = obs.sd)

y3.test.impl = impl(em = y3.test.pred$mean, em.sd = y3.test.pred$sd,
                    disc = disc, obs = y3.test, disc.sd = disc.sd, obs.sd = obs.sd)


# Now, how similar are the NROY spaces?
# Pairs plot with the truth
dfunc.up <- function(x,y,...){
  # function for plotting 2d kernel density estimates in pairs() plot.
  require(MASS)
  require(RColorBrewer)
  
  br <- brewer.pal(9, 'Blues')
  kde <- kde2d(x,y)
  image(kde, col = br, add = TRUE)
}

dfunc.up.truth = function(x,y, ...){
  # function for plotting 2d kernel density estimates in pairs() plot,
  # adding a data point overlay.
  require(MASS)
  require(RColorBrewer)
  
  xtrue <- tail(x,1)
  ytrue <- tail(y,1)
  
  xdash <- head(x, -1)
  ydash <- head(y, -1)
  
  br <- brewer.pal(9, 'Blues')
  
  kde <- kde2d(xdash,ydash)
  image(kde, col = br, add = TRUE)
  points(xtrue, ytrue, pch =21, col = 'black', bg = 'red', cex = 1.5)
}

pairs(rbind(X.unif[bt.ix, ], x.test), panel = dfunc.up.truth, gap = 0, upper.panel = NULL)


#intersections
# These are different lengths. That makes things tricky because of course,
# different outputs behave differently due to each input. 
# This might scupper the whole idea.

y1.bt.ix = which(y1.impl<thres)
y2.bt.ix = which(y2.impl<thres)
y3.bt.ix = which(y3.impl<thres)

X.nroy.y1 = X.unif[y1.bt.ix, ]
X.nroy.y2 = X.unif[y2.bt.ix, ]
X.nroy.y3 = X.unif[y3.bt.ix, ]

test = rbind(X.nroy.y1, X.nroy.y2, X.nroy.y3)

cols = c(rep('black', length(y1.bt.ix)), rep('red', length(y2.bt.ix)),rep('blue', length(y3.bt.ix)))

pairs(test, col = cols)


sorensen_dice = function(a,b){
  
  out = (2 * (length(intersect(a,b)))) / (length(a) + length(b))
  out
  
}

sorensen_dice(y1.bt.ix,y2.bt.ix )
sorensen_dice(y2.bt.ix,y3.bt.ix )
sorensen_dice(y1.bt.ix,y3.bt.ix )

# Different outputs suggest different parts of parameter space, and that's normal.
# Part of the power of this process in the augmented paper is that we'd expect the
# three observations to suggest similar parts of parameter space.

# We still haven't got far in deciding if something is a disrepancy, or just that 
# a different part of parameter space would work better.

plot(X.unif[,4], y1.impl)


# TO DO:
# Multi-wave history matching.
# Adding a discrepancy and finding it.


# create a mmlh within the bound of the (original?) input space
# See which members have max(I) < threshold and keep them.
# Run the model at those points and rebuild the emulator
# Sample uniformly across original space and see what proportion
# of inputs are below the threshold.






fit.list = list(fit.y1, fit.y2, fit.y3)
obs.list = list(y1.test, y2.test, y3.test)
obs.sd.list = list(0,0,0)
disc.list = list(0,0,0)
disc.sd.list = list(0.3, 0.3, 0.3)


generate.nroy = function(fit.list, n.cand =10000, obs.list, thres = 3){
  # Generate design samples that are in nroy space.
  
  # Set up some holders
  impl.mat = matrix(NA, nrow = n.cand, ncol = length(fit.list))
  
  X.cand = maximinLHS(n = n.cand, k = k, dup = 2)
  X.cand = maximinLHS(n = n.cand, k = k, dup = 2)
  
  # just loop it to start with
  for(i in 1:length(fit.list)){
    
    fit = fit.list[[i]]
    obs = obs.list[[i]]
    obs.sd = obs.sd.list[[i]]
    disc = disc.list[[i]]
    disc.sd = disc.sd.list[[i]]
    
    pred = predict(fit, newdata = X.cand, type = 'UK')
    
    pred.impl = impl(em = pred$mean, em.sd = pred$sd,
                     disc = disc, obs = obs, disc.sd = disc.sd, obs.sd = obs.sd)
    
    impl.mat[,i]  = pred.impl
  }
  
  # helper function
  all.bt = function(x, thres) all(x < thres)
  # where is below the threshold for all cases?
  nroy.ix = apply(impl.mat, 1, FUN = all.bt, thres = thres)
  
  X.nroy = X.cand[nroy.ix, ]
  X.nroy
}


# this is very slow if you use mmlh and a large set of candidate points (say 10000)
nroy.test = generate.nroy(fit.list, n.cand = 10000, obs.list = obs.list, thres = 3)

Y.wave2 = matrix(NA, nrow = nrow(nroy.test), ncol = 3)

for(i in 1:nrow(nroy.test)){
  
  x = nroy.test[i, ]
  
  y1.wave2 = friedman(x)
  y2.wave2 = park1(x)
  y3.wave2 = park2(x)
  
  Y.wave2[i, ] = c(y1.wave2, y2.wave2, y3.wave2)
}

# These are the augmented designs - with new design points from the initially NROY
# space.
X.aug2 = rbind(X, nroy.test)
Y.aug2 = rbind(Y, Y.wave2)


fit2.y1 = km(~., design = X.aug2, response = Y.aug2[, 1])
fit2.y2 = km(~., design = X.aug2, response = Y.aug2[, 2])
fit2.y3 = km(~., design = X.aug2, response = Y.aug2[, 3])

# sample uniformly from input space
#n.unif = 1000
#X.unif = samp.unif(n.unif, mins = rep(0,k), maxes = rep(1,k))
#colnames(X.unif) <- c('x1', 'x2', 'x3', 'x4')

# predict the model output at each sample point
y1.pred2 = predict(fit2.y1, newdata = X.unif, type = 'UK')
y2.pred2 = predict(fit2.y2, newdata = X.unif, type = 'UK')
y3.pred2 = predict(fit2.y3, newdata = X.unif, type = 'UK')

y1.impl2 = impl(em = y1.pred2$mean, em.sd = y1.pred2$sd,
               disc = disc, obs = y1.test, disc.sd = disc.sd, obs.sd = obs.sd)

y2.impl2 = impl(em = y2.pred2$mean, em.sd = y2.pred2$sd,
               disc = disc, obs = y2.test, disc.sd = disc.sd, obs.sd = obs.sd)

y3.impl2 = impl(em = y3.pred2$mean, em.sd = y3.pred2$sd,
               disc = disc, obs = y3.test, disc.sd = disc.sd, obs.sd = obs.sd)

# candidate mmlh point generation

# CAN AUGMENT THE LHS - DO THIS!
# Could generate an initial large mmlh and just subsample (hmm)
# Augment each time in waves



# Find the implausibility of a large sample from the emulator.

# sample uniformly from input space
n.unif = 1000
X.unif = samp.unif(n.unif, mins = rep(0,k), maxes = rep(1,k))
colnames(X.unif) <- c('x1', 'x2', 'x3', 'x4')

# predict the model output at each sample point for each fit in the list

y1.pred = predict(fit.y1, newdata = X.unif, type = 'UK')
y2.pred = predict(fit.y2, newdata = X.unif, type = 'UK')
y3.pred = predict(fit.y3, newdata = X.unif, type = 'UK')

unif.pred.list = function(fit.list, X.unif){
  
  
}







