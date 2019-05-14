%% Set up problem parameters.
k = 4;
%p = [500, 1000];

%p = 25*2.^(0:8);
p = [3000:500:5000, 6000:1000:10000];
%p = [150:50:500, 600:100:1000, 1250: 250:2500];
%p = 2500;

p = 1500;

N = 25*ones(k,length(p));
blocksize = ceil(p/(2*k));
Ntest = 500*ones(k, length(p));
T = 1;
savemat = true;
r = 0.5;

%% Run experiment.

[times, errs, feats]=time_compare_1(p,r,k,blocksize, N,Ntest, T, savemat);
fprintf('DONEZO!\n\n')


%% Save.
save 'TIMEEXPTSRESULTS.mat'

%% SCRATCH.
