

% Basic Variable Declarations

t = [-12:.01:12]; % plot mesh
p = .5; % sampling rate
N = 10000; % approximation cutoff
L = 100; % erasure set size
W = 100; % bridge set size

% La = [1:1:L]+N+1; % erasure set
% Om = [-(L/2-1):0,L+1:3/2*L]+N+1; % bridge set

% La = [0:1:1];
% La = [La,La+4,La+8,La+12,La+16,La+20,La+24,La+28,La+32,La+36,La+40,La+44,La+48,La+52,La+56,La+60,La+64,La+68,La+72,La+76]+N+1;
% Om = [2:1:3];
% Om = [Om,Om+4,Om+8,Om+12,Om+16,Om+20,Om+24,Om+28,Om+32,Om+36,Om+40,Om+44,Om+48,Om+52,Om+56,Om+60,Om+64,Om+68,Om+72,Om+76]+N+1;

% La = [1:5:5*L]+N+1;
% Om = [0:5:5*L-1,2:5:5*L+1]+N+1;

La = [1:2:2*L-1]+N+1;
Om = [2:2:2*L]+N+1;

% Om = [1:1:L]+N+1; % erasure set
% La = [-(L/2-1):0,L+1:3/2*L]+N+1; % bridge set

% Computing the sampling coefficients

FC = zeros(1,2*N+1); % sampling coefficients
for(n = -N:1:N)
  FC(n+N+1) = sinc(pi * n * p); % sampling coefficients of sinc(pi*x)
end
FC1 = FC; % actual sampling coefficients
FC(La) = zeros(size(La)); % sampling coefficients with erasures
FC2 = FC;

% Computing the N-term approximation of f_R

f_Rapprox = zeros(size(t)); % N-term approximation of f_R
for(n = -N:1:N)
  f_Rapprox = f_Rapprox + FC(n+N+1) * p * sinc(pi * (t-n*p));
end

% Computing the actual f_R

f_R = sinc(pi*t); % actual f_R
for(n = La)
  f_R = f_R - FC1(n) * p * sinc(pi * (t-(n-N-1)*p));
end

% Computing the Bridge Matrix

B = zeros(L,W); % bridge matrix
for(j = 1:1:L)
  for(k = 1:1:W)
    B(j,k) = sinc(pi*p*(Om(k)-La(j)));
  end
end

% condy = cond(B)

% Computing the Right hand side of the bridge equation

RHS = zeros(L,L); % right hand side
for(j = 1:1:L)
  for(k = 1:1:L)
    RHS(j,k) = sinc(pi*p*(La(k)-La(j)));
  end
end

% Solving for C

C = B \ RHS; % coefficient matrix

norm(C)

% Reconstructing the sampling coefficients

CfROapprox = zeros(size(Om)); % sampling coefficients of f_R over the bridge set.
for(n=-N:1:N)
  CfROapprox = CfROapprox + p*FC(n+N+1)*sinc(pi*(p*(Om-N-1)-p*n));
end

CfRLapprox = zeros(size(La)); % sampling coefficients of f_R over the erasure set.
for(n=-N:1:N)
  CfRLapprox = CfRLapprox + p*FC(n+N+1)*sinc(pi*(p*(La-N-1)-p*n));
end

CfRO = sinc(pi*p*(Om-N-1));
for(n=La)
  CfRO = CfRO - p*FC1(n)*sinc(pi*(p*(Om-N-1)-p*(n-N-1)));
end

CfRL = sinc(pi*p*(La-N-1));
for(n=La)
  CfRL = CfRL - p*FC1(n)*sinc(pi*(p*(La-N-1)-p*(n-N-1)));
end

FC(La) = (C' * (FC(Om)' - CfROapprox') + CfRLapprox')'; % implementation of the reconstruction algorithm.

FC2(La) = (C' * (FC2(Om)' - CfRO') + CfRL')';

N_term_approx_error = max(abs(FC1-FC)) % accuracy check
error = max(abs(FC1-FC2))

f_Bapprox = zeros(size(t)); % N term approximation of reconstructed function
for(n = -N:1:N)
  f_Bapprox = f_Bapprox + p* FC(n+N+1)*sinc(pi*(t-n*p));
end

f_B = f_R; % the reconstructed function
for(n = La-N-1)
  f_B = f_B + p* FC2(n+N+1)*sinc(pi*(t-n*p));
end

f_E = sinc(pi*t)-f_R;
f_E_no_approx = f_B - f_R;
f_E_approx = f_Bapprox - f_Rapprox;

% plot(t,f_E,'k');
% hold on;
% plot(t,f_E_no_approx,'g');
% plot(t,f_E_approx,'r');
% hold off;

figure;
plot(t,f_Bapprox,'r');
hold on;
xlim([-12,12]);
ylim([-1.1,1.1]);
%plot(t,f_Rapprox,'r');
%plot(t,f_R,'g')
plot(t,sinc(pi*t),'--b');
% plot(t,f_B,'g');
legend('Bridged Function','Original Function','Location','southeast')
hold off;