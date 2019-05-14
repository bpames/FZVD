%% Set up problem parameters.
k = 4;
%p = [500, 1000];

%p = 25*2.^(0:8);
%p = [50:50:500 , 600:100:1000, 1250: 250:2500, 3000:500:4000];

%p = [150:50:500, 600:100:1000, 1250: 250:2500];
p = 2500;
N = 25*ones(k,length(p));
blocksize = ceil(p/(2*k));
Ntest = 500*ones(k, length(p));
T = 5;
savemat = true;
r = 0.5;

%% Run experiment.



[times, errs, feats]=time_compare_1(p,r,k,blocksize, N,Ntest, T, savemat);
fprintf('DONEZO!\n\n')
% %% Plots.
% figure; hold on
% set(gca,'XScale','log','YScale','log');
% for i = 1:4
%     loglog(p, mean(times(:,:,i)))
% end
% legend('ball', 'sphere', 'asda', 'classic')
% hold off
% 
% %% Error plots.
% figure; hold on
% % set(gca,'XScale','log','YScale','log');
% for i = 1:4
%     plot(p, mean(errs(:,:,i)))
% end
% legend('ball', 'sphere', 'asda', 'classic')
% hold off
% 
% %% Feature plots.
% figure; hold on
% %set(gca,'XScale','log','YScale','log');
% for i = 1:4
%     plot(p, mean(feats(:,:,i)))
% end
% legend('ball', 'sphere', 'asda', 'classic')
% hold off



%% Save.
save 'TIMEEXPTSRESULTS.mat'

%% SCRATCH.
