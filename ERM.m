% Parameters:
% m = Height of Erasure Recovery Matrix.
%     m must be larger than size of erasure set.
% n = Dimension
% N = Length of the Frame
% L = The Erasure Set

m = 250;
n = 250;
N = 1000;
L = [1:10];

% This block creates a Dual Frame pair (DF,EF), with erasure recovery matrix M
% according to the GEGW construction algorithm for erasure recovery matrices.

A = randn(N,2*n+m);
[A,~] = qr(A,0);
DF = sqrt(N/n)*A(:,1:n)';
EF = sqrt(n/N)*A(:,n+1:2*n)' + (n/N)*DF;
M = sqrt(N/m)*A(:,2*n+1:2*n+m)';

% f is a random vector that we will try to recover 
% from frame coefficient erasures.

f = randn(n,1);
f = f./norm(f);

% FC are the frame coefficients of f.
        
FC = EF' * f;

% We erase the frame coefficients indexed by
% L.

FC(L) = zeros(size(L'));

% We compute f_R.

f_R = DF*FC;

% We reconstruct the erased frame coefficients.
        
LC = setdiff(1:N,L);
FC(L) = -(M(:,L)' * M(:,L))\(M(:,L)' * (M(:,LC) * FC(LC)));

% We reconstruct the signal.

g = f_R + DF(:,L) * FC(L);

% We compute the \ell^2 norm of the reconstruction
% error.
        
norm(f-g)
