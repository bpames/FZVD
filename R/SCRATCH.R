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

cmns <- penzda(Xt = Xtrain, Yt = Ytrain, type ="sphere")

