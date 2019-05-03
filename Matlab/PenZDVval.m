function [val_w, DVs, gamma,gammas, its, w, scaler, val_score, classMeans] = PenZDVval(train, val,D,num_gammas,gmults, sparsity_level, beta, tol, maxits,quiet)
%Atrain: training data
%D: Penalty dictionary basis matrix
%Aval: validation set
%k: #of classes within the training and validation sets
%num_gammas: number of gammas to train on
%g_mults:(c_min, c_max): parameters defining range of gammas to train g_max*(c_min, c_max)

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

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%Initialize the validation scores
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
val_score=zeros(num_gammas,1);
%mc_ind=1;
%l0_ind=1;
best_ind=1;
%min_mc=1;
%min_l0=p;
triv=0;

%Initialize DVs and iterations
N0=N;
DVs=zeros(p,K-1,num_gammas);
its=zeros(K-1,num_gammas);

%For each gamma, calculate ZVDs and corresponding validation score
gammas=zeros(num_gammas,K-1);

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Cross-validation.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Set initial gamma.
gammas(:,1)=gamscale*norm(RN*w,2)^2/norm(Dx(N*w),1);

for i=1:num_gammas
    N=N0;
    %R=R0;
    
    % Set gamma_1.
    gammas(i,1) = gmult(i)*norm(RN*w,2)^2/norm(Dx(N*w,1)); 
    
    for j=1:K-1
        
        % Initial solutions.
        if i == 1
            sols0.x = w;
            sols0.y = Dx(N*sols0.x);
            sols0.z = zeros(p,1);
        else
            sols0.y = DVs(:,j,i-1); % Warm-start with previous DV.
            sols0.x = N'*Dtx(sols0.y); % x = N'*D'*y.
            sols0.z = zeros(p,1);
        end
        quietADMM=1;
        
        %Call ADMM
        
        if isequal(consttype,'ball')
            
            % Call ball-constrained solver.
            [~,y,~,its]=SZVD_ADMM_V2(R,N,RN, D,sols0,gammas(i,j),beta,tol,maxits,quietADMM);
            
            % Normalize y, if necessary.
            DVs(:,j,i) = y/norm(y);
            
        elseif isequal(consttype,'sphere')
            
            % Call spherical solver.
            [~,DVs(:,j,i),~,tmpits]=SZVD_ADMM_S(R,N,RN, D,sols0,gammas(i,j),beta,tol,maxits,quietADMM);
            
        else % ERROR.
            error('Invalid constraint type. Please indicate if using inequality ("ball") or equality ("sphere") constraints.')
        end
          
        its(j,i)=tmpits;
        %Update N and B 
        if j< K-1
            %Project N onto orthogonal complement of Nx 
            x=Dx(N*x);
            x=x/norm(x);
            N=Nupdate1(N,x);
            RN = R*N;
            [~,sigma,w]=svds(RN, 1, 'largest');
            R=R/sigma;
            % Set gamma.
            gammas(i,j+1)=gmult(i)*norm(RN*w,2)^2/norm((D*N*w),1);
        end
    end
    %Get performance scores on the validation set
    %Call test ZVD to get predictions
    [stats,~,proj,cent]=predict(DVs(:,:,i), val, classMeans);%%%%%%%%%%%%%
    proj;
    cent;
    %If gamma gives trivial sol, give it a large penalty
    if (sum(stats.l0)<3)
        val_score(i)=100*size(val,1);
        triv=1;
    else
        if (sum(stats.l0)>sparsity_level*size(DVs(:,:,i),1))
            val_score(i)=sum(stats.l0);
        else
            val_score(i)=stats.mc;
        end
    end
    %Update the best gamma so far
    if (val_score(i)<=val_score(best_ind))
        best_ind=i;
    end
    %Record sparest nontrivial sol 
    %if (min(stats.l0)>3 && stats.l0<min_l0)
    %    l0_ind=1;
    %    l0_x=DVs(:,:,i);
    %    min_10=stats.l0;
    %end
    %Record best misclassificartion error
    %if(stats.mc<=min_mc)
    %    min_ind=i;
    %    mc_x=DVs(:,:,i);
    %    min_mc=stats.mc;
    %end
    %Terminate if a trivial solution has been found
    if (quiet==0)
       fprintf('it = %g, val_score= %g, mc=%g, l0=%g, its=%g \n', i, val_scores(i), stats.mc, sum(stats.l0), mean(its(:,i)));
    end
    if(triv==1)
        break;
    end
end
%Export DVs found using validation
val_w=DVs(:,:,best_ind);
gamma=gammas(best_ind,:);
scaler=gmults(best_ind);
val_score=val_score(best_ind);
