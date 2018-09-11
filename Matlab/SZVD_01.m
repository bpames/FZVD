function [DVs,its,pen_scal,N,classMeans,mus]=SZVD_01(train,gamma,D,penalty,scaling,tol,maxits,beta,quiet)
%Prepare the training data (ZVD step)
get_DVs=1;
[~,p]=size(train); 
p=p-1;
w0 = ZVD(train, scaling, get_DVs);
classMeans=w0.means;
mus=w0.mu;
if (penalty==1)
    s=sqrt(diag(w0.W));
else
    s=ones(p,1);
end
w0.s=s;
%w0.B = 1/2*(w0.B + w0.B');
%B0 = w0.N' * w0.B * w0.N;
d1 = w0.dvs(:,1);
B0 = (w0.N' * w0.B * w0.N)/(d1'*w0.B*d1);
B0 = (B0+B0')/2;
N = w0.N;
K=w0.k;
%Initialization for the output
DVs=zeros(p,K-1);
its=zeros(1,K-1);

% Define d operators.
if isdiag(D)
    Dx = @(x) diag(D).*x; % diagonal scaling is D is diagonal.
    Dtx = @(x) diag(D).*x; 
else
    Dx = @(x) D*x;
    Dtx = @(x) D'*x;
end

%Call ADMM

% Initial solutions.
sols0.x = N'*Dtx(w0.dvs(:,1));
sols0.y = w0.dvs(:,1);
sols0.z = zeros(p,1);

for i=1:(K-1)
    %Initial solutions.
   
    [x,~,~,its]=SZVD_ADMM(B0,N,D,sols0,s,gamma(i),beta,tol,maxits,quiet);
    
    DVs(:,i)=Dx(N*x);
    DVs(:,i) = DVs(:,i)/norm(DVs(:,i));
    its(i)=its;
    if (quiet == 0)          
        fprintf('Found SZVD %g after %g its \n', i, its(i));
    end
   
    if(i<(K-1))
        %Project N onto orthogonal complement of Nx 
        x=DVs(:,i);
        x=x/norm(x);
        Ntmp=N-x*(x'*N);%null complement of W
        %QR factorization of Ntmp
        [Q,R]=qr(Ntmp);
        %Extract non-zero rows of R to get columns of Q to update N
        R_rows=(abs(diag(R))>1e-6);
        N=Q(:,R_rows);
        B0=N'*w0.B*N;
        B0=0.5*(B0+B0');
        
        
        % Update initial solutions.
        [w,~] = eigs(B0, 1, 'LM');
        
        %plot(w)
        % Project back to the original space.
        d1 = N * w;
        
        sols0.y = Dx(d1);
        sols0.x = w;        
        sols0.z = zeros(p,1);
        
        % Rescale B0.
        B0 = B0/(d1'*w0.B*d1);
        
        %eigs(B0)
    end
   
end
pen_scal=s;
end

        
        
    
    



