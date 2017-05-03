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

M = randn(N,m);
[M,~] = qr(M,0);
M = sqrt(N/m) * M';
A = [M',randn(N,n)];
[A,~] = qr(A,0);
F = A(:,m+1:m+n)';

% f is a random vector that we will try to recover 
% from frame coefficient erasures.

f = randn(n,1);
f = f./norm(f);

% FC are the frame coefficients of f.
        
FC = F' * f;

% We erase the frame coefficients indexed by
% L.

FC(L) = zeros(size(L'));

% We compute f_R.

f_R = F*FC;

% We reconstruct the erased frame coefficients.
        
LC = setdiff(1:N,L);
FC(L) = -(M(:,L)' * M(:,L))\(M(:,L)' * (M(:,LC) * FC(LC)));

% We reconstruct the signal.

g = f_R + F(:,L) * FC(L);

% We compute the \ell^2 norm of the reconstruction
% error.
        
norm(f-g)