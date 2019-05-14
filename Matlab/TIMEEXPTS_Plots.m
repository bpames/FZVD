%% Plots.
figure; hold on
set(gca,'XScale','log','YScale','log');
for i = 1:4
    loglog(p, mean(times(:,:,i)))
end
legend('ball', 'sphere', 'asda', 'classic')
hold off

%% Error plots.
figure; hold on
% set(gca,'XScale','log','YScale','log');
for i = 1:4
    plot(p, mean(errs(:,:,i)))
end
legend('ball', 'sphere', 'asda', 'classic')
%xlim([300, 2500])
herold off

%% Feature plots.
figure; hold on
%set(gca,'XScale','log','YScale','log');
for i = 1:4
    plot(p, mean(feats(:,:,i)))
end
legend('ball', 'sphere', 'asda', 'classic')
hold off

