clc;

% Parameters:
% n = Dimension
% N = Length of the Frame
% L = The Erasure Set

n = 2000;
N = 3000;
L = [1:900];

% The columns of F are a randomly generated frame.
% The columns of G are the standard dual to F.

F = rand(n,N);
S = F * F';
G = S \ F;

% f is a random vector that we will try to recover 
% from frame coefficient erasures.

f = rand(n,1);
f = f ./ norm(f,2);

% FC are the frame coefficients of f.

FC = G' * f;

% We erase the frame coefficients indexed by
% L, 

FC(L) = zeros(size(L'));

% We compute f_R.

f_R = F * FC;

% We compute the matrix M and the coefficient
% matrix.

M = (F(:,L)' * G(:,L))';
C = (eye(max(size(L))) - M) \ eye(max(size(L)));

% We compute the reconstruction.

H = G(:,L)' * f_R;
g = f_R;
g = g + F(:,L) * C * H;

% We compute the \ell^2 norm of the error in the
% reconstruction.

norm(f-g,2)