% SCRATCH.

%% Load data.
clc
clear
% load('OOdata (normalized).mat')
%load('Coffee-normalized.mat')
 load('ECGdata (normalized).mat')
[n,p] = size(train);
p = p-1;


%% Split testing data as validation and testing set.
valratio = 0.35; % Extract this fraction for validation set.
[val, test] = train_test_split(test, valratio);



%%
%prepare the data set
gmults = linspace(0, 1.5, 15);
beta=2;
tol.rel = 1e-3;
tol.abs= 1e-3;
maxits= 1000;
quiet=false;

consttype = 'ball';
consttype = 'sphere';

sparsity_level = 0.25;

D = eye(p);

%% Call validation set.
[val_w, DVs, gamma, best_ind, val_score, classMeans] = PenZDAval(train, val,D, gmults, consttype, sparsity_level, beta, tol, maxits,quiet);
gamma, best_ind

stats = predict(val_w, test, classMeans)

% Plot DVs.
[~,K] = size(classMeans);
for i = 1:K-1
    figure
    plot(1:p, val_w(:,i))
end

%% Try cross-validation.

nfolds = 10;
[bestDVs, bestind, inds, bestgamma,  cv_scores, classMeans] = PenZDAcv(train, nfolds, D, gmults, consttype, sparsity_level, beta, tol, maxits,quiet);

bestind
bestgamma
cv_scores
cvstats = predict(bestDVs, test, classMeans)

% Plot DVs.
[~,K] = size(classMeans);
for i = 1:K-1
    figure
    plot(1:p, bestDVs(:,i))
end

%% Call ASDA.

% Make Y.
[nt, p] = size(train);
Xt = train(:, 2:p);
p = p-1;

labs = train(:,1);
Yt = zeros(nt, max(labs));
for i = 1:nt
    Yt(i, labs(i)) = 1;
end

%% 
Om = eye(p);
gam = 0.1;
lam = 0.15; % What is a better choice?
q = 1;
PGsteps = 500;
PGtol = 1e-3;
maxits = 500;
tol = 1e-3;

[B,Q] = SDAAP(Xt, Yt, Om, gam, lam, q, PGsteps, PGtol, maxits, tol);
ASDAstats = predict(B, test, classMeans)


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