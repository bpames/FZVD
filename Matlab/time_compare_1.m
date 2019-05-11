function [times, errs, feats]=time_compare_1(p,r,k,blocksize, N,Ntest, T, savemat)
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
scaling = 1;
penalty=0;
beta=3;
tol.rel = 1e-3;
tol.abs= 1e-3;
maxits=500;
quiet=1;


%Initialize matrices for storing results
times = zeros(T, length(p),4);
errs = times;
feats = times;

% Set up timing table
fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
fprintf('p   \t | Ball \t | Sphere \t | ASDA \t | Old \n')
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
        
        % Reset method number.
        meth = 0;
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % SOLVE USING BALL CONSTRAINED PENZDA.
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        meth = meth + 1;
        consttype = 'ball';        
        tic;
        [DVs,~, ~,~,classMeans, gamma] = PenZDA(train,D, tol,maxits,beta,quiet, consttype,gamscale);
              
        times(j, i, meth) =toc; % Stop timer after training is finished.        
        
        % Sample and normalize test data.
        test = type1_data(p(i),r,k,Ntest(:, i), blocksize(i)); 
        test_obs=test(:,2:(p(i)+1));
        test_obs=normalize_test(test_obs,mu_train,sig_train);
        test(:, 2:(p(i)+1)) = test_obs ;
        
        % Check classification and feature selection performance.
        [stats,~,~,~]=predict(DVs,test,classMeans);
        errs(j,i, meth)=stats.mc;
        feats(j,i, meth)=sum(stats.l0);
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % SOLVE USING SPHERICALLY CONSTRAINED PROBLEM.
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        meth = meth + 1;
        consttype = 'sphere';        
        tic;
        [DVs,~, ~,~,classMeans,~] = PenZDA(train,D, tol,maxits,beta,quiet, consttype,gamscale);              
        times(j,i, meth) =toc; % Stop timer after training is finished.        
        % Check classification and feature selection performance.
        [stats,~,~,~]=predict(DVs,test,classMeans);
        errs(j,i, meth)=stats.mc;
        feats(j,i, meth)=sum(stats.l0);
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % SOLVE USING SDAD or SDAAP.
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        meth = meth + 1;
        
        % Make Y.
        nt = size(train,1);
        Xt = train(:, 2:(p(i) + 1));        
        
        labs = train(:,1);
        Yt = zeros(nt, k);
        for ii = 1:nt
            Yt(ii, labs(ii)) = 1;
        end
        
        % Set ASDA parameters.
        tmp = rand(p(i));
        Om = D + 1e-2*(tmp*tmp');
        gam = 0.1;
        %lam = 0.15; % What is a better choice?
        q = k-1;
        PGsteps = 500;
        PGtol = 1e-3;
        maxits = 500;
        ASDAtol = 1e-3;
        
        % Call ASDA-APG.
        tic
        % Add classMeans calculation.
        classes=train(:,1);
        [nt,pt]=size(train);
        X=train(:,2:pt);
        %X=normalize(X);
        %Extract observations
        labels=unique(classes);
        K=length(labels);
        %Initiate matrix of within-class means
        pt=pt-1;
        classMeans=zeros(pt,K);
        ClassMeans=zeros(pt,K);
        
        % Initialize gamma.
        gamma = zeros(K-1,1);
        
        %for each class, make an object in the list containing only the obs of that
        %class and update the between and within-class sample
        M=zeros(pt,nt);
        
        for ii=1:K
            class_obs=X(classes==labels(ii),:);
            %Get the number of obs in that class (the number of rows)
            ni=size(class_obs,1);
            %Compute within-class mean
            classMeans(:,ii)=mean(class_obs);
            %Update W
            xj=class_obs-ones(ni,1)*classMeans(:,ii)';
            M(:,classes == labels(ii)) =xj';
            ClassMeans(:,ii)=mean(class_obs)*sqrt(ni);
        end
    
        % Add lambda calculation.
        A = 2*(Xt'*Xt + gam*Om);
        % Precompute Mj = I - Qj*Qj'*D.
        Qj = ones(K, 1);
        Di = 1/nt*(Yt'*Yt);
        Mj = @(u) u - Qj*(Qj'*(Di*u));
        
        % Initialize theta.
        theta = Mj(rand(K,1));
        theta = theta/sqrt(theta'*Di*theta);
        
        %%
        % Form d.
        d = 2*Xt'*Yt*theta/nt;
        
        % Initialize beta.
        beta = A\d; % 1st unpenalized solution.
        
        % Choose lambda so that unpenalized solution always has negative value.
        lmax = (beta'*d - 0.5*beta'*A*beta)/norm(beta, 1);
        
        % Set lambda.
        lam = gamscale*lmax;

        % Call SDAAP.
        [DVs,~] = SDAAP(Xt, Yt, Om, gam, lam, q, PGsteps, PGtol, maxits, ASDAtol);
        times(j,i,meth) = toc;
        [stats,~,~,~]=predict(DVs,test,classMeans);
        errs(j,i, meth)=stats.mc;
        feats(j,i, meth)=sum(stats.l0);
        
        
        % Calculate gamma/lambda.
        
        
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        % SOLVE USING OLD CODE.
        %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        meth = meth + 1;
        beta = 3;
        %Repeat using the old code and save results to remaining matrices.
        tic;
        [DVs2,~,~,~,~,~]=SZVD_01(train,gamma,D,penalty,scaling,tol,maxits,beta,1);

        times(j,i, meth) =toc;
        
        %fprintf('old-test')
        [stats, ~] = predict(DVs2,test,classMeans);
        errs(j,i, meth)=stats.mc;
        feats(j,i, meth)=sum(stats.l0);
        
        % Print intermediate stats.
        fprintf('%4d \t | %1.3f \t | %1.3f \t | %1.3f \t | %1.3f \n', p(i), times(j,i, 1), times(j,i, 2), times(j,i, 3), times(j,i, 4))
        
        % Save workspace if desired.
        if(savemat)
            save 'timecompareres.mat' times errs feats
        end
    end
end





