% SCRATCH.

%% Load data.
clc
clear
load('OOdata (normalized).mat')
[n,p] = size(train);
p = p-1;

%%
%prepare the data set
gamscale=1.5;
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
pentype = 'ball';
[DVs,~,~,~,~,classMeans, ~] = PenZDA(train,D,penalty,tol,maxits,beta,quiet, pentype,gamscale);
             
t0 = toc, % Stop timer after training is finished.
        
stats0 = test_ZVD_V1(DVs,test,classMeans)

for i = 1:3
    figure
    plot(DVs(:,i))
end


%% Call spherical solver.
tic
gamscale = 1;
pentype = 'sphere';
[DVs1,~,~,~,~,classMeans,gamma] = PenZDA(train,D,penalty,tol,maxits,beta,quiet, pentype,gamscale);
t1 = toc, % Stop timer after training is finished.
        
stats1 = test_ZVD_V1(DVs1,test,classMeans)

for i = 1:3
    figure
    plot(DVs1(:,i))
end


%% Check accuracy.

t0 - t1
norm(DVs - DVs1)