function [val_w, DVs, gamma, best_ind, val_score, classMeans] = PenZDAval(train, val,D, gmults, consttype, sparsity_level, beta, tol, maxits,quiet)
% PENZDAVAL performs validation to train regularization parameters and DVs.
%
% INPUT.
% train, val: training and validation set.
% D: Penalty dictionary basis matrix 
% gmults: list of multipliers of maximum gamma (calculated later) to compare.
% consttype: "ball/sphere" imposes ball/spherical constraint on DVs.
% sparsity_level: desired minimum sparsity level for validation scoring.
% beta: augmented Lagrangian penalty parameter.
% tol: stopping tolerances
% maxits: maximum number of iterations.
% quiet = true suppresses display of intermediate outputs.
% OUTPUT.
% val_w: discriminant vectors corresponding to best parameter.
% DVs: set of all discriminant vectors.
% gamma: best parameter choice.
% best_ind: position of best gamma.
% val_score: array of validation scores.
% classMeans: set of training data class means.

% Initialize classMeans of training data and get number of classes.
% classes=train(:,1);
[n,p]=size(train);
p = p-1;
% X=train(:,2:p);
% %X=normalize(X);
% %Extract observations 
% labels=unique(classes); 
% K=length(labels);
% %Initiate matrix of within-class means
% p=p-1;
% classMeans=zeros(p,K);
% R=zeros(K,p);
% 
% %for each class, make an object in the list containing only the obs of that
% %class and update the between and within-class sample
% M=zeros(p,n);
% 
% 
% % Calculate scaled classMeans for calculating R and N in objective function.
% for i=1:K    
%     class_obs=X(classes==labels(i),:);
%     %Get the number of obs in that class (the number of rows)
%     ni=size(class_obs,1);
%     %Compute within-class mean
%     classMeans(:,i)=mean(class_obs);
%     %Update W 
%     xj=class_obs-ones(ni,1)*classMeans(:,i)';
%     M(:,classes == labels(i)) =xj';
%     R(i,:)= sqrt(ni)*mean(class_obs)';
% end

% Call calcClassMeans to compute class-means and covariance matrices.
[classMeans, K, M, R]=calcClassMeans(train);

% Calculate null basis.
N=null(M);

%Compute leading eigenvector of N'*B*N
RN = R*N;
%size(RN)
[~,sigma,w]=svds(RN, 1,'largest');
% normalize R.
R=R/sigma;
RN = RN/sigma;
RN0 = RN;
Nw = N*w;

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
num_gammas = length(gmults);
val_score= (p+1)*ones(num_gammas,1);
%mc_ind=1;
%l0_ind=1;
best_ind=1;
%min_mc=1;
%min_l0=p;
triv=0;

%Initialize DVs and iterations
N0=N;
DVs=zeros(p,K-1,num_gammas);
its=zeros(num_gammas,1);

%For each gamma, calculate ZVDs and corresponding validation score
gammas=zeros(num_gammas, K-1);
gmax = norm(RN*w,2)^2/norm(Dx(Nw), 1); 


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Validation step.
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Set initial gamma.
%gammas(:,1)=gamscale*norm(RN*w,2)^2/norm(Dx(N*w),1);

for i=1:num_gammas
    N=N0;
    RN = RN0;
    
    % Set gamma_1.
%     size(RN)
%     size(N)
%     size(w)    
    gammas(:, 1) = gmults(i)*gmax;    
    
    for j=1:K-1
        
        % Initial solutions.
        if i == 1
            sols0.x = w;
            sols0.y = Dx(Nw);
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
            [x,y,~,tmpits]=SZVD_ADMM_V2(R,N,RN, D,sols0,gammas(i,j),beta,tol,maxits,quietADMM);
            
            % Normalize y, if necessary.
            DVs(:,j,i) = y/norm(y);
            
        elseif isequal(consttype,'sphere')
            
            % Call spherical solver.
            [x,DVs(:,j,i),~,tmpits]=SZVD_ADMM_S(R,N,RN, D,sols0,gammas(i,j),beta,tol,maxits,quietADMM);
            
        else % ERROR.
            error('Invalid constraint type. Please indicate if using inequality ("ball") or equality ("sphere") constraints.')
        end
          
        its(i)= its(i) + tmpits;
        %Update N and B 
        if j< K-1
            %Project N onto orthogonal complement of Nx 
            x=Dx(N*x);
            x=x/norm(x);
%             fprintf('update N\n')
            N=Nupdate1(N,x);
            RN = R*N;
            [~,sigma,w]=svds(RN, 1, 'largest');
            R=R/sigma;
            % Set gamma.
%             size(RN)
%             size(N)
%             size(w)
            gammas(i,j+1)=gmults(i)*norm(RN*w,2)^2/norm(Dx(N*w),1);
        end
    end
    %Get performance scores on the validation set
    %Call test ZVD to get predictions
    [stats,~,~,~]=predict(DVs(:,:,i), val, classMeans);%%%%%%%%%%%%%
    
    %If gamma gives trivial sol, give it a large penalty
    if (min(stats.l0)<1)
        val_score(i)=(p+1)*size(val,1);
        triv=1;
    else
        if (sum(stats.l0)>sparsity_level*(K-1)*p)
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
    
    if (quiet==false)
       if (i==1) 
           fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
           fprintf('i \t\t + gam \t\t\t + score \t\t + mc \t\t\t + l0 \t\t + its \n')
           fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
       end
       fprintf('%3d \t +  %1.2e \t + %1.2e \t + %1.2e \t + %1.2e \t + %1.2e \n', i, gammas(i), val_score(i),  stats.mc, sum(stats.l0), its(i));
    end
    if(triv==1)
        break;
    end
end


%Export DVs found using validation
val_w=DVs(:,:,best_ind);
gamma=gammas(best_ind,:);
%scaler=gmults(best_ind);
%val_score=val_score(best_ind);




