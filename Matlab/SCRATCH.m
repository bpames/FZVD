% SCRATCH.

%% Load data.
clc
clear
% load('OOdata (normalized).mat')
% load('Coffee-normalized.mat')
load('ECGdata (normalized).mat')
[n,p] = size(train);
p = p-1;


%% Split testing data as validation and testing set.
valratio = 0.05; % Extract this fraction for validation set.
[val, test] = train_test_split(test, valratio);



%%
%prepare the data set
gmults = linspace(0, 1.5, 15);
beta=2;
tol.rel = 1e-5;
tol.abs= 1e-5;
maxits=100;
quiet=false;

consttype = 'sphere';

sparsity_level = 0.3;

D = eye(p);

%% Call validation set.
[val_w, DVs, gamma,gammas, its, w, scaler, val_score, classMeans] = PenZDAval(train, val,D, gmults, consttype, sparsity_level, beta, tol, maxits,quiet);
gamma

stats = predict(val_w, test, classMeans)

%% Call solver.
gamscale = 0.3;
tic;
pentype = 'ball';
[DVs,~,~,~,classMeans, ~] = PenZDA(train,D,tol,maxits,beta,quiet, pentype,gamscale);
             
t0 = toc, % Stop timer after training is finished.
        
stats0 = test_ZVD_V1(DVs,test,classMeans)

% [~,K] = size(classMeans);
% for i = 1:K-1
%     figure
%     plot(DVs(:,i))
% end
% 

% Call spherical solver.
tic
pentype = 'sphere';
[DVs1,~,~,~,classMeans,gamma] = PenZDA(train,D, tol,maxits,beta,quiet, pentype,gamscale);
t1 = toc, % Stop timer after training is finished.
        
stats1 = test_ZVD_V1(DVs1,test,classMeans)

% Plot DVs.
[~,K] = size(classMeans);
for i = 1:K-1
    figure
    plot(1:p,DVs1(:,i),  1:p,DVs(:,i) )
end


%% Check accuracy.

t0 - t1
norm(DVs - DVs1)