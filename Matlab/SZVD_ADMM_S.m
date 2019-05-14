function [x,y,z,its, errtol] = SZVD_ADMM_S(R, N, RN, D, sols0,gamma, beta, tol, maxits, quiet)

% Iteratively solves the problem
%       min{-1/2*x'B'x + gamma p(y): l2(y) = 1, DNx = y}
% using ADMM.
%====================================================================
% Input.
%====================================================================
%   R: Matrix stores the class means, i.e. R=classMeans'.
%   N: basis matrix for null space of covariance matrix W.
%   RN: product of RN used repeatedly.
%   D: penalty dictionary/basis.
%   sols0: initial solutions sols0.x, sols0.y, sols0.z
%   pen_scal: penalty scaling term.
%   gamma:  l1 regularization parameter.
%   beta:    penalty term controlling the splitting constraint.
%   tol:    tol.abs = absolute error, tol. rel = relative error to be
%                   achieved to declare convergence of the algorithm.
%   maxits: maximum number of iterations of the algorithm to run.
%   quiet: toggles between displaying intermediate statistics.
%====================================================================
% Output:.
%====================================================================
%   x, y, z: iterates at termination.
%   its: number of iterations required to converge.
%   errtol: stopping error bound at termination.


%====================================================================
% Precomputes quantities that will be used repeatedly by the algorithm.
%====================================================================

% Dimension of decision variables.
p = size(D, 1);

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

% RN product.
% RN = R*N;
% size(RN)

%====================================================================
% Initialize x solution and constants for x update.
%====================================================================
% K > 2 Case.
K= size(R, 1);
% Compute initial x.
x = sols0.x;
Nx = N*x;
 % Take Cholesky of beta I - B (for use in update of x)
%V=chol(eye(K)-1/beta*R*(N*N')*R','upper');
%[V1,V2] = qr(eye(K)-1/beta*R*(N*N')*R');
[P,L] = lu(eye(K)-1/beta*(RN*RN'));


%====================================================================
% Initialize decision variables y and z.
%====================================================================

% Initialize y and z.
y = sols0.y;
z = sols0.z;

%====================================================================
%% Call the algorithm.
%====================================================================

for iter=1:maxits    
    
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % Step 1: Perform shrinkage to update y_k+1.
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    % Save previous iterate.
    yold = y;
    
    
    % Calculate largest magnitude entry of b.
    b = Dx(Nx) + z;
    [mx, ix] = max(abs(b));
    
    % Update y. 
    if mx <= gamma  % all-zeros is optimal for ball-constrained problem.
        % Set yix to be indicator for largest value.
        y = zeros(p,1);
        y(ix) = sign(b(ix));
    else % all-zeros is suboptimal. Use soft-threshholding.
        y = vec_shrink(b, gamma);
        y = y/norm(y); % Normalize y (if necessary).
    end
    
    % Truncate complex part (if have +0i terms appearing.)
    y = real(y);
    
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % Step 2: Update x_k+1 by solving
    % x_k+1 = argmin { -x'*A*x + beta/2 l2(x - y_k+1 + z_k)^2}
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    % Compute RHS.
    b = N'*Dtx(beta*y - z);
    
    % K > 2.
    % Update using by solving the system LL'x = b.
    xtmp=P\(RN*b);
    xtmp=L\(xtmp);
    x=1/beta*b+1/beta^2*RN'*xtmp;
    
    % Truncate complex part (if have +0i terms appearing.)
    x = real(x);
    
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %  Step 3: Update the Lagrange multipliers
    % (according to the formula z_k+1 = z_k + beta*(N*x_k+1 - y_k+1) ).
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %zold = z;
    % Primal residual.
    Nx = N*x;
    r = Dx(Nx) - y;    
    
    % Update z.
    z = real(z + beta*r);
    
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % Check stopping criteria.
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    %----------------------------------------------------------------
    % Primal constraint violation.
    %----------------------------------------------------------------
    
    % l2 norm of the residual.
    dr = norm(r);
	
	%----------------------------------------------------------------   
    % Dual constraint violation.
    %----------------------------------------------------------------
    % Dual residual.
    s = beta*(y - yold);    
    % l2 norm of the residual.
    ds = norm(s);
    
	%----------------------------------------------------------------
    % Check if the stopping criteria are satisfied.
	%----------------------------------------------------------------
    
    % Compute absolute and relative tolerances.
    ep = sqrt(p)*tol.abs + tol.rel*max(norm(x), norm(y));
    es = sqrt(p)*tol.abs + tol.rel*norm(y);

    
    % Display current iteration stats.
    if (quiet==0 && mod(iter, 1) == 0)
        fprintf('it = %g, primal_viol = %3.2e, dual_viol = %3.2e, norm_DV = %3.2e\n', iter, dr-ep, ds-es, norm(y))
    end
    
    % Check if the residual norms are less than the given tolerance.
    if (dr < ep && ds < es && iter > 10)
        break % The algorithm has converged.
    end
    
end %for.


%====================================================================
% Output results.
%====================================================================

if maxits > 0
    its = iter;
    errtol = min(ep,es);
else
    its = 0;
    errtol=0;
end

end
