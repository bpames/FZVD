function [DVs,x, its,pen_scal,N,classMeans, gamma]=PenZDA(train,D,penalty,tol,maxits,beta,quiet,gamscale)

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

%Compute leading eigenvector of N'*B*N
RN = R*N;
%size(RN)
[~,sigma,w]=svds(RN, 1,'largest');
% normalize R.
R=R/sigma;
RN = RN/sigma;

% Define d operators.
if isdiag(D)
    Dx = @(x) diag(D).*x; % diagonal scaling is D is diagonal.
    Dtx = @(x) diag(D).*x; 
else
    Dx = @(x) D*x;
    Dtx = @(x) D'*x;
end


% Set gamma.
gamma(1)=gamscale*norm(RN*w,2)^2/norm(Dx(N*w),1);

% ppt = toc;
% fprintf('ppt %1.4d \n', ppt)

%Set sigma--penalty parameter
if (penalty==1)
    s=sqrt(diag(M*M'));
else
    s=ones(1,p);
end
%Initialization for the output
DVs=zeros(p,K-1);
its=zeros(1,K-1);
%Call ADMM

for i=1:(K-1)
    %Initial solutions.
%     tic
    sols0.x = w;
    sols0.y = Dx(N*w);
    sols0.z = zeros(p,1);
    [x,y,~,its]=SZVD_ADMM_V2(R,N,RN, D,sols0,s,gamma(i),beta,tol,maxits,quiet);
    DVs(:,i) = y/norm(y);
%     st = st + toc;
%     fprintf('solve time %1.4d \n', st)
    %DVs(:,i)=D*N*x;
    its(i)=its;
    if (quiet == 0)          
        fprintf('Found SZVD %g after %g its \n', i, its(i));
    end
    
%     tic
    if(i<(K-1))
        %Project N onto orthogonal complement of Nx 
        x=Dx(N*x);
        x=x/norm(x);
        N=Nupdate1(N,x);
        RN = R*N;
        [~,sigma,w]=svds(RN, 1, 'largest');
        R=R/sigma;
        % Set gamma.
        gamma(i+1)=gamscale*norm(RN*w,2)^2/norm((D*N*w),1);
    end
%     ntime = ntime + toc;
%     fprintf('Nt %1.4d \n', ntime)
    
end
pen_scal=s;
end