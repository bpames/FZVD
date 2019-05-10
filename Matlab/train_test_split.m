function [train,test]=train_test_split(traindata,trainratio)
N=size(traindata,1);
tf = false(N,1);  
tf(1:round(trainratio*N)) = true;     
tf = tf(randperm(N));   
train= traindata(tf,:); 
test = traindata(~tf,:); 