function [bestDVs, bestind, bestgamma,  cv_scores, classMeans] = PenZDAcv(train, nfolds, D, gmults, consttype, sparsity_level, beta, tol, maxits,quiet)
%PENZDACV cross-validation for choice of parameters in PenZDA.
%
% INPUT.
% train: training data.
% nfolds: number of folds for cross-validation.
% D: Penalty dictionary basis matrix 
% gmults: list of multipliers of maximum gamma (calculated later) to compare.
% consttype: "ball/sphere" imposes ball/spherical constraint on DVs.
% sparsity_level: desired minimum sparsity level for validation scoring.
% beta: augmented Lagrangian penalty parameter.
% tol: stopping tolerances
% maxits: maximum number of iterations.
% quiet = true suppresses display of intermediate outputs.
%
% OUTPUT.
% bestDVs: discriminant vectors corresponding to best parameter.
% bestind: position of multiplier for best gamma.
% bestgamma: best parameter choice.
% cv_scores: array of cross-validation scores.

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Initialize CV scores and intermediate outputs.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% number of gammas to test for each fold.
ngam = length(gmults);

% Array of cross-validation scores.
cv_scores = zeros(ngam, nfolds);

% number and dimension of training observations.
[n,p] = size(train); p = p - 1;

% Set training ratio. (Ignored if nfolds = n for LOO).)
tratio = 1- nfolds/n;

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% CV scheme.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if (quiet == false)
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    fprintf('Leave-One-Out CV\n')
end

for f = 1:nfolds
    
    %=============================================================
    % Split data.
    %=============================================================
    
    if (nfolds == n) % Leave-one-out CV.        
        % Leave out f-th observation as validation observation.        
        trainf = train; trainf(f,:) = [];
        valf = train(f,:);
%         size(trainf)
%         size(valf)
    elseif (nfolds > n) % error.
        error('number of folds cannot exceed number of training observations.')
    else % use nfolds/n observations for validation.              
       % Split training set into training data and validation data.
       [trainf, valf] = train_test_split(train, tratio);              
    end
    
    % Normalize training data
    xf = trainf(:, 2:(p+1));
    [xnf, muf, sigf] = normalize(xf);
    trainf(:, 2:(p+1)) = xnf;
    
    % "Normalize" validation data.
    valf(:,2:(p+1)) = normalize_test(valf(:,2:(p+1)),muf,sigf);

    %=============================================================
    % Calculate validation scores for this fold.
    %=============================================================
    % Display fold number.
    if (quiet == false)
        fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
        fprintf('Fold %3d\n', f)        
    end
    
    % Call PenZDAval.
    [~, ~, ~, ~, cv_scores(:, f), ~] = PenZDAval(trainf, valf,D, gmults, consttype, sparsity_level, beta, tol, maxits,quiet);    


end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Choose "best" multiplier to have best average score.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Calculate average classification scores.
meanscores = mean(cv_scores,2); %size(meanscores)

% Find position of best score.
[bestscore,~] = min(meanscores);
inds = find(meanscores == bestscore);
bestind = max(inds);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Solve for discriminant vectors using the full training set and best
% parameters.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
gamscale = gmults(bestind);
[bestDVs,~, ~,~,classMeans, bestgamma] = PenZDA(train,D, tol,maxits,beta,1, consttype, gamscale);

