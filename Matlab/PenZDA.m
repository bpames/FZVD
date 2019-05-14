function [DVs,x, its,N,classMeans, gamma]=PenZDA(train,D, tol,maxits,beta,quiet, consttype,gamscale)

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
R=zeros(K,p);

% Initialize gamma.
gamma = zeros(K-1,1);

%for each class, make an object in the list containing only the obs of that
%class and update the between and within-class sample
M=zeros(n,p);

for i=1:K    
    class_obs=X(classes==labels(i),:);
    %Get the number of obs in that class (the number of rows)
    ni=size(class_obs,1);
    %Compute within-class mean
    classMeans(:,i)=mean(class_obs);
    %Update W and R.
    M(classes == labels(i),:) =class_obs-ones(ni,1)*classMeans(:,i)';
    R(i,:)= sqrt(ni)*classMeans(:,i)';
end

%Find null basis. 
N=null(M);

%Compute leading eigenvector of N'*B*N
RN = R*N;
[~,sigma,w]=svds(RN, 1,'largest');
Nw = N*w;

% normalize R.
R=R/sigma;
RN = RN/sigma;

% Define d operators.
if isdiag(D)
    % Check if D = I
    d = diag(D); 
    if norm(d - ones(p,1)) < 1e-12 % D=I
        Dx = @(x) x;
        %Dtx = Dx;
    else % D ~=I
        Dx = @(x) d.*x; % diagonal scaling if D is diagonal.
        %Dtx = Dx;
    end
else
    Dx = @(x) D*x;
%     Dtx = @(x) D'*x;
end


% Set gamma.
gamma(1)=gamscale*norm(RN*w,2)^2/norm(Dx(Nw),1);

%Initialization for the output
DVs=zeros(p,K-1);
its=zeros(1,K-1);

%Call ADMM
for i=1:(K-1)
    %Initial solutions.
%     tic
    sols0.x = w;
    sols0.y = Dx(Nw);
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
         
        % Calculate discriminant vector from x.
%         x=Dx(N*x);     
% x=x/norm(x);
        
        % Try y instead.
        x = DVs(:,i);
        
        
        
        % Call N update.
        N=Nupdate1(N,x);       
        
        % Update RN and initial solution.
        RN = R*N;
        [~,sigma,w]=svds(RN, 1, 'largest');
        
        R=R/sigma;
        % Set gamma.
        Nw = N*w;
        gamma(i+1)=gamscale*norm(RN*w,2)^2/norm(Dx(Nw),1);
    end
    
end

end