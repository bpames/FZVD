function [obs,mu,sigma]=type1_data(p,r,k,N, blocksize)
%p: number of features; p=k*l
%r: value of constant covariance between features
%k: number of classes
%N: vector of number of observations per class
%blocksize: size of each block of values distinguishing class-means.

%set size of diagonal blocks of sigma


%initialize mean vectors mu and covariance matrix sigma
mu=zeros(p,k);
 

%mu_ij=0.7 if blocksize(i-1)+ 1 <=j<=blocksize*i
for i=1:k
    mu(((i-1)*blocksize+1):(blocksize*i),i)= ones(blocksize,1);
    %mu(1:l,i)=0.7*ones(l,1);
end


%form the covariance matrix sigma
sigma=r*ones(p,p);
sdiag = diag(sigma);
sigma = eye(p) + sigma - diag(sdiag);
%sigma(logical(eye(p)))=1;

%generate sum(N) observations using mu and sigma
%initialize observation matrix
obs=zeros(sum(N),(p+1));

%initialize the start/end of current class
start=1;
ending=N(1);

%generate i-th class
for i=1:k
    %sample observations
    obs(start:ending,1)=i*ones(N(i),1);
    %size(mu)
    obs(start:ending, 2:(p+1))=mvnrnd(mu(:,i),sigma,N(i));
    
    %update block positions
    if i<k
        start=ending+1;
        ending=ending+N(i+1);
    end
end







