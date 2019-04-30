function N1=Nupdate1(N,v)
% NUPDATE1 Updates null basis using single Householder reflection.
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT.
% N - current null-basis matrix.
% v - new discriminant vector.
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT.
% N1 - matrix with columns giving orthonormal basis for intersection of
% col(N) and orthogonal complement of v.
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Extract column of partial QR factorization (used to find N) to reflect
% onto x-axis.
x = N'*v;

% Extract dimension of current null-space.
d = length(x);

% Update x to (implicitly) calculate Householder reflector.
x(1) = sign(x(1))*norm(x) + x(1);
x = x./norm(x);

% Update N using x.
N1 = N(:, 2:d) - 2*(N*x)*x(2:d)';