# Make ASDA Input.
n <- dim(train)[1]
p <- dim(train)[2] - 1

Xtrain <- as.matrix(train[, 2:(p+1)])
Ytrain <- as.factor(train$X1)

Xtest <- as.matrix(test[, 2:(p+1)])
Ytest <- as.factor(test$X1)

ntest <- dim(test)[1]



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Call ASDA.
library(accSDA)

res <- ASDA(Xt = Xtrain, Yt = Ytrain)

# Predict on the test data
preds <- predict(res, newdata = Xtest)

# Calculate accuracy
sum(preds$class == Ytest)/ntest # We have N samples per class, so total 3*N

# Plot DV.
plot(res$beta, type = 'l')


#++++++++++++++++++++++++++++++++++++++++++++++++++++
# Test penzda.
library(MASS)
library(rARPACK)

cmns <- penzda(Xt = Xtrain, Yt = Ytrain, maxits=50, tol = 1e-3, type ="ball")
cmns$DVs
plot(cmns$DVs[,1], type="l")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# TEST VALIDATION SCHEME.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Split data.
set.seed(2)
nval <- round(0.2*n)
valinds <-  sample.int(n, size = nval, replace = FALSE)

Xval <- Xtrain[valinds, ]
Yval <- Ytrain[valinds]
Xvt <- Xtrain[-valinds, ]
Yvt <- Ytrain[-valinds]

# Need to renormalize training and validation data.
trainlist <- normalize(x = Xvt)
Xvt <- trainlist$x
mu <- trainlist$mu
sig <- trainlist$sig

Xval <- normalizetest(x=Xval, mu = mu, sig = sig)

res <- penzdaVAL(Xt = Xvt, Yt = Yvt, Xval = Xval, Yval=Yval, maxits = 500,
                 gmults = c(0.25, 0.5, 0.75, 1),sparsity_level = 0.4,
                 quiet = FALSE,type = "sphere")
                    

plot(res$DVs, type="l")
res$val_score

penstats <- predict(obj = res, Xtest = Xtest, Ytest = Ytest)
penstats$mc

# TEST VALIDATION SCHEME.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

cvres <- penzdaCV(Xt = Xtrain, Yt = Ytrain, nfolds = 5, maxits = 500,
                 gmults = c(0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2,5),sparsity_level = 0.4,
                 quiet = FALSE,type = "sphere")

plot(cvres$DVs, type="l")
cvres$cvscores

rowMeans(cvres$cvscores)

penstats <- predict(obj = cvres, Xtest = Xtest, Ytest = Ytest)

penstats$mc
penstats$l0/136

cbind(penstats$preds, Ytest, penstats$dist)

# v  <- c(1, 2, -4)
# a <- 1.5
# 
# s <- vecshrink(v,a)

#++++++++++++++++++++++++++++++++
x <- cmns$sols$x
y <- cmns$sols$y
z <- cmns$sols$z
gamma <- cmns$gamma
D <- diag(p)



# r <- rbind(c(1,2,3),
#            c(0,2,1),
#            c(0,0,2))
# 
# prds <- max.col(r)
# r

# y <- backsolve(r, x <- c(8,4,2)) 
# 
# r %*% y
# 
# backsolve(r, x, transpose = TRUE)

