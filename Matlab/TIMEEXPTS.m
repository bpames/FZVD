%% Set up problem parameters.
k = 4;
%p = [500, 1000];

%p = 25*2.^(0:8);
<<<<<<< HEAD
p = [50:50:500 , 600:100:1000, 1250: 250:2500, 3000:500:4000];
=======
p = [50:50:500];%, 600:100:1000, 1250: 250:2500, 3000:500:5000];
>>>>>>> parent of 8db2e06... Optimize SZVD_V6
N = 10*ones(k,length(p));
blocksize = ceil(p/(2*k));
Ntest = 500*ones(k, length(p));
T = 100;
savemat = true;
r = 0.5;

%% Run experiment.
[time1,time2, err1, err2, feat1, feat2]=time_compare_1(p,r,k,blocksize, N,Ntest, T, savemat);

<<<<<<< HEAD
%% Plots.
%loglog(1:length(p), mean(time1), 1:length(p), mean(time2))

save('TIMEEXPTSRESULTS.mat')
=======
%% SCRATCH.
i = 2; j =1;

% Generate  and process (i,j)th training data.
train = type1_data(p(i),r,k,N(:, i), blocksize(i));
[train_obs, mu_train, sig_train] = normalize(train(:,2:(p(i)+1)));
train=[train(:,1), train_obs];
        
 % Solve using new version and record cpu time.
        size(train)
        tic;
        [DVs,~,~,~,~,classMeans,gamma] = SZVD_V6(train,D,penalty,tol,maxits,beta,quiet,gammascale);
              
        time1(j, i) =toc; % Stop timer after training is finished.
        

%%
p = p(i);

%%
% tic
classes=train(:,1);
[n,p]=size(train);
X=train(:,2:p);
%X=normalize(X);
%Extract observations 
labels=unique(classes); 
K=length(labels);
%Initiate matrix of within-class means
p=p-1;
classMeans=zeros(p,K);
ClassMeans=zeros(p,K);

% Initialize gamma.
gamma = zeros(K-1,1);

%for each class, make an object in the list containing only the obs of that
%class and update the between and within-class sample
M=zeros(p,n);

for i=1:K    
    class_obs=X(classes==labels(i),:);
    %Get the number of obs in that class (the number of rows)
    ni=size(class_obs,1);
    %Compute within-class mean
    classMeans(:,i)=mean(class_obs);
    %Update W 
    xj=class_obs-ones(ni,1)*classMeans(:,i)';
    M(:,classes == labels(i)) =xj';
    ClassMeans(:,i)=mean(class_obs)*sqrt(ni);
end

%Symmetrize W and B
R=ClassMeans';
%Find ZVDs 
N=null(M');

    %Compute k-1 nontrivial e-vectors of N'*B*N
    RN = R*N;
    %size(RN)
    [~,sigma,w]=svd(RN);
    size(RN)
    w=w(:,1);
    % normalize R.
    R=R/sigma(1,1);
    RN = RN/sigma(1,1);
    
    % Set gamma.
    gamma(1)=gamscale*norm(RN*w,2)^2/norm((D*N*w),1);
>>>>>>> parent of 8db2e06... Optimize SZVD_V6

