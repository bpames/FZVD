function [bestDVs, allDVs,] = PenZDAcv(train, D, gmults, consttype, sparsity_level, beta, tol, maxits,quiet)
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
% val_w: discriminant vectors corresponding to best parameter.
% bestgamma: best parameter choice.
% bestind: position of best gamma.
% cv_scores: array of cross-validation scores.

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Initialize CV scores and intermediate outputs.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ngam = length(gmults);
cv_scores = zeros(ngam, nfolds);



end

