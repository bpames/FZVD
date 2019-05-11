%% Set up problem parameters.
k = 4;
%p = [500, 1000];

%p = 25*2.^(0:8);
%p = [50:50:500 , 600:100:1000, 1250: 250:2500, 3000:500:4000];

p = [50:50:250];%, 600:100:1000, 1250: 250:2500, 3000:500:5000];
N = 10*ones(k,length(p));
blocksize = ceil(p/(2*k));
Ntest = 500*ones(k, length(p));
T = 3;
savemat = true;
r = 0.5;

%% Run experiment.
[times, errs, feats]=time_compare_1(p,r,k,blocksize, N,Ntest, T, savemat);

%% Plots.
figure; hold on
set(gca,'XScale','log','YScale','log');
for i = 1:4
    loglog(p, mean(times(:,:,i)))
end
hold off

%% Save.
save 'TIMEEXPTSRESULTS.mat'
%% SCRATCH.
