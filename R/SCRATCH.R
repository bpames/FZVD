# Make ASDA Input.
n <- dim(train)[1]
p <- dim(train)[2] - 1

Xtrain <- as.matrix(train[, 2:(p+1)])
Ytrain <- train$X1

Xtest <- as.matrix(test[, 2:(p+1)])
Ytest <- test$X1

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

cmns <- penzda(Xt = Xtrain, Yt = Ytrain, maxits=15, tol = 1e-3, type ="ball")
cmns$DVs
plot(cmns$DVs, type="l")

penstats <- predict(obj = cmns, Xtest = Xtest, Ytest = Ytest)

v  <- c(1, 2, -4)
a <- 1.5

s <- vecshrink(v,a)

#++++++++++++++++++++++++++++++++
x <- cmns$sols$x
y <- cmns$sols$y
z <- cmns$sols$z
gamma <- cmns$gamma
D <- diag(p)



r <- rbind(c(1,2,3),
           c(0,2,1),
           c(0,0,2))

prds <- max.col(r)
r

y <- backsolve(r, x <- c(8,4,2)) 

r %*% y

backsolve(r, x, transpose = TRUE)

