% SCRATCH.

%% Load data.
clc
clear
load('OOdata (normalized).mat')
[n,p] = size(train);
p = p-1;

%%
%prepare the data set
gammascale=0.15;
penalty=0;
scaling=0;
beta=2;
tol.rel = 1e-5;
tol.abs= 1e-5;
maxits=100;
quiet=0;

D = eye(p);

%% Call solver.
tic;
[DVs,~,~,~,~,classMeans,gamma] = SZVD_V6(train,D,penalty,tol,maxits,beta,quiet,gammascale);
             
t0 = toc, % Stop timer after training is finished.
        
stats0 = test_ZVD_V1(DVs,test,classMeans)

plot(DVs)


%% Call new solver.