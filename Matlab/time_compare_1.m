function [timeB,timeS, time0,  ErrB, ErrS, Err0, FeatB, FeatS, Feat0]=time_compare_1(p,r,k,blocksize, N,Ntest, T, savemat)
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
gamscale=0.5;
penalty=0;
scaling=1;
beta=3;
tol.rel = 1e-3;
tol.abs= 1e-3;
maxits=100;
quiet=1;


%Initialize matrices for storing results
timeB = zeros(T, length(p));
time0=timeB;ErrB=timeB;Err0=timeB;FeatB=timeB;Feat0 = timeB;
timeS=timeB; ErrS = timeB; FeatS = timeB;

% Set up timing table
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
fprintf('p   \t | Ball \t | Sphere \t | Old \n')
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')

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
        
   
        % SOLVE USING BALL CONSTRAINED PENZDA.
        consttype = 'ball';        
        tic;
        [DVs,~, ~,~,classMeans, gamma] = PenZDA(train,D, tol,maxits,beta,quiet, consttype,gamscale);
              
        timeB(j, i) =toc; % Stop timer after training is finished.        
        
        % Sample and normalize test data.
        test = type1_data(p(i),r,k,Ntest(:, i), blocksize(i)); 
        test_obs=test(:,2:(p(i)+1));
        test_obs=normalize_test(test_obs,mu_train,sig_train);
        test(:, 2:(p(i)+1)) = test_obs ;
        
        % Check classification and feature selection performance.
        [stats,~,~,~]=predict(DVs,test,classMeans);
        ErrB(j,i)=stats.mc;
        FeatB(j,i)=sum(stats.l0);
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % SOLVE USING SPHERICALLY CONSTRAINED PROBLEM.
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        consttype = 'sphere';        
        tic;
        [DVs,~, ~,~,classMeans, gamma] = PenZDA(train,D, tol,maxits,beta,quiet, consttype,gamscale);              
        timeS(j, i) =toc; % Stop timer after training is finished.        
        % Check classification and feature selection performance.
        [stats,~,~,~]=predict(DVs,test,classMeans);
        ErrS(j,i)=stats.mc;
        FeatS(j,i)=sum(stats.l0);
        
        
                   
        %Repeat using the old code and save results to remaining matrices.
        tic;
        [DVs2,~,~,~,~,~]=SZVD_01(train,gamma,D,penalty,scaling,tol,maxits,beta,1);

        time0(j, i) =toc;
        
        %fprintf('old-test')
        [stats, ~] = predict(DVs2,test,classMeans);
        Err0(j,i)=stats.mc;
        Feat0(j,i)=sum(stats.l0);
        
        % Print intermediate stats.
        fprintf('%4d \t | %1.4f \t | %1.4f \t | %1.4f \n', p(i), timeB(j,i), timeS(j,i), time0(j,i))
        
        % Save workspace if desired.
        if(savemat)
            save 'timecompareres.mat' timeB timeS time0 ErrB ErrS Err0 FeatB FeatS Feat0
        end
    end
end





