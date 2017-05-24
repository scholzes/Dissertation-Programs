
% Basic Variable Declarations

t = [-12:.01:12]; % plot mesh
p = .5; % sampling rate
N = 10000; % approximation cutoff
L = 12; % erasure set size

La = [1:1:L]+N+1; % erasure set
% La = [1:2:2*L-1]+N+1;

% Computing the sampling coefficients

FC = zeros(1,2*N+1); % sampling coefficients
for(n = -N:1:N)
  FC(n+N+1) = sinc(pi * n * p); % sampling coefficients of sinc(pi*x)
end
FC1 = FC; % actual sampling coefficients
FC(La) = zeros(size(La)); % sampling coefficients with erasures

% Computing the N-term approximation of f_R

f_R = zeros(size(t)); % N-term approximation of f_R
for(n = -N:1:N)
  f_R = f_R + FC(n+N+1) * p * sinc(pi * (t-n*p));
end

% Computing the Bridge Matrix

M = zeros(L,L); % bridge matrix
for(j = 1:1:L)
  for(k = 1:1:L)
    M(j,k) = p * sinc(pi*p*(La(k)-La(j)));
  end
end

% Solving for C

C = (eye(L) - M) \ (eye(L)); % coefficient matrix

norm(C)

% Reconstruction

CfRL = zeros(size(La)); % sampling coefficients of f_R over the erasure set.
for(n=-N:1:N)
  CfRL = CfRL + p*FC(n+N+1)*sinc(pi*(p*(La-N-1)-p*n));
end

f_B = f_R; % the reconstructed function
for(j = 1:1:L)
    for(k = 1:1:L)
        f_B = f_B + p * C(j,k) * CfRL(k) * sinc(pi*(t-p*(La(j)-N-1)));
    end
end

figure;
plot(t,f_B,'r');
hold on;
xlim([-12,12]);
ylim([-0.75,1.1]);
%plot(t,f_Rapprox,'r');
%plot(t,f_R,'g')
plot(t,sinc(pi*t),'--b');
% plot(t,f_B,'g');
legend('Reconstructed Function','Original Function','Location','northeast')
hold off;