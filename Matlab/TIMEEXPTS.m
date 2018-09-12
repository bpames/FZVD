%% Set up problem parameters.
k = 4;
%p = [500, 1000];

%p = 25*2.^(0:8);
%p = [50:50:500 , 600:100:1000, 1250: 250:2500, 3000:500:4000];

p = [50:50:500];%, 600:100:1000, 1250: 250:2500, 3000:500:5000];
N = 10*ones(k,length(p));
blocksize = ceil(p/(2*k));
Ntest = 500*ones(k, length(p));
T = 10;
savemat = true;
r = 0.5;

%% Run experiment.
[time1,time2, err1, err2, feat1, feat2]=time_compare_1(p,r,k,blocksize, N,Ntest, T, savemat);


%% Plots.
loglog(1:length(p), mean(time1), 1:length(p), mean(time2))

save 'TIMEEXPTSRESULTS.mat'
%% SCRATCH.
