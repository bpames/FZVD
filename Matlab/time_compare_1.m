function [time1,time2, NumErr1, NumErr2, NumFeat1, NumFeat2]=time_compare_1(p,r,k,blocksize, N,Ntest, T, savemat)
%p: vector of number of features
%r: value of constant covariance between features
%k: number of classes
%blocksize: size of each block of values distinguishing class-means.
%N: array of vectors of number of training observations per class.
%Ntest: array of vector of number of testing observations per class.
%T: number of trials for p
%savemat: logical indicating whether to save intermediate workspace to
%file.

%prepare the data set
gammascale=0.5;
penalty=0;
scaling=1;
beta=3;
tol.rel = 1e-3;
tol.abs= 1e-3;
maxits=100;
quiet=1;


%Initialize matrices for storing results
time1 = zeros(T, length(p));
time2=time1;NumErr1=time1;NumErr2=time1;NumFeat1=time1;NumFeat2 = time1;

% Set up timing table
fprintf('+++++++++++++++++++++++++++++++++++\n')
fprintf('p   | New     | Old \n')
fprintf('+++++++++++++++++++++++++++++++++++\n')

for i=1:length(p)
    % Update dictionary matrix.
    D=eye(p(i));
    
    %+++++++++++++++++++++++++++++++++++++++
    % Run trials for the ith problem size.
    for j=1:T
        
        % Generate  and process (i,j)th training data.                        
        train = type1_data(p(i),r,k,N(:, i), blocksize(i)); 
        [train_obs, mu_train, sig_train] = normalize(train(:,2:(p(i)+1)));
        train=[train(:,1), train_obs];
        
%         train=[train(:,1),normalize(train(:,2:end))];
%         train_obs=train(:,2:end);
%         mu_train=mean(train_obs);
%         sig_train=std(train_obs);
        
    %fprintf('new')
    
        % Solve using new version and record cpu time.
        %size(train)
        tic;
        [DVs,~,~,~,~,classMeans,gamma] = SZVD_V6(train,D,penalty,tol,maxits,beta,quiet,gammascale);
              
        time1(j, i) =toc; % Stop timer after training is finished.
        
        %fprintf('new-test')
        % Sample and normalize test data.
        test = type1_data(p(i),r,k,Ntest(:, i), blocksize(i)); 
        test_obs=test(:,2:(p(i)+1));
        test_obs=normalize_test(test_obs,mu_train,sig_train);
        test(:, 2:(p(i)+1)) = test_obs ;
        
        % Check classification and feature selection performance.
        [stats,~,~,~]=test_ZVD_V1(DVs,test,classMeans);
        NumErr1(j,i)=stats.mc;
        NumFeat1(j,i)=sum(stats.l0);
                   
        %Repeat using the old code and save results to remaining matrices.
        tic;
        [DVs2,~,~,~,~,~]=SZVD_01(train,gamma,D,penalty,scaling,tol,maxits,beta,1);

        time2(j, i) =toc;
        
        %fprintf('old-test')
        [stats, ~] = test_ZVD_V1(DVs2,test,classMeans);
        NumErr2(j,i)=stats.mc;
        NumFeat2(j,i)=sum(stats.l0);
        
        % Print intermediate stats.
        fprintf('%4d | %1.4f | %1.4f \n', p(i), time1(j,i), time2(j,i))
        
        % Save workspace if desired.
        if(savemat)
            save 'timecompareres.mat' time1 time2 NumErr1 NumErr2 NumFeat1 NumFeat2
        end
    end
end





