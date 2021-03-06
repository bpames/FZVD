function [DVs,x, its,N,classMeans, gamma]=GENZDA(train,D, tol,maxits,beta,quiet, Q, pentype, consttype, alpha, gamscale)

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
    % Check if D = I
    d = diag(D); 
    if norm(d - ones(p,1)) < 1e-12 % D=I
        Dx = @(x) x;
        Dtx = Dx;
    else % D ~=I
        Dx = @(x) d.*x; % diagonal scaling if D is diagonal.
        Dtx = Dx;
    end
else
    Dx = @(x) D*x;
    Dtx = @(x) D'*x;
end


% Set gamma.
gamma(1)=gamscale*norm(RN*w,2)^2/norm(Dx(N*w),1);

% ppt = toc;
% fprintf('ppt %1.4d \n', ppt)

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
    
    if isequal(consttype,'ball')
        
        % Call ball-constrained solver.
        [x,y,~,its]=SZVD_ADMM_V2(R,N,RN, D,sols0,gamma(i),beta,tol,maxits,quiet);
        
        % Normalize y, if necessary.
        DVs(:,i) = y/norm(y);
    
    elseif isequal(consttype,'sphere')
        
        % Call spherical solver.
        [x,DVs(:,i),~,its]=SZVD_ADMM_S(R,N,RN, D,sols0,gamma(i),beta,tol,maxits,quiet);
        
    else % ERROR. 
        error('Invalid constraint type. Please indicate if using inequality ("ball") or equality ("sphere") constraints.')
    end
    
    
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

end