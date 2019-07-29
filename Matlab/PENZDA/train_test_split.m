function [train,test]=train_test_split(data,trainratio)
% TRAIN_TEST_SPLIT divides data into training and validation sets.
% INPUT.
% data - data set to split.
% trainratio - fraction to be included in training set.
% OUTPUT.
% train, test - data split into training and validation sets.


% Check for valid inputs.
if (trainratio > 1 || trainratio < 0)
    error('trainratio must be in [0,1]')
end

% Get number of observations.
N=size(data,1);

% Initialize flags for training (T) and validation (F) sets.
tf = false(N,1);  
tf(1:round(trainratio*N)) = true;     

% Randomly permute tf to get list of training observations.
tf = tf(randperm(N));   

% Split data according to tf.
train= data(tf,:); 
test = data(~tf,:); 