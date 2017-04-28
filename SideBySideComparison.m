clc;

% Original Images are 256 pixels X 256 pixels.

fprintf('Reading Image \n');

COMPRESSION_PERCENT = 0.03; % Compressed Signal will be approximately
% n = 256^2 * COMPRESSION_PERCENT dimensional.
snr = .05;
percenterasures = .05;

Original_Image_Double = double(imread('Lena.bmp'));

fprintf('Performing Image Compression \n');

Compressed_Image_Double = fft(reshape(Original_Image_Double,[256*256,1]));
[S,I] = sort(abs(Compressed_Image_Double),'descend');
n = round(COMPRESSION_PERCENT*256*256);
LSC = Compressed_Image_Double(I(n+1:256*256));
Compressed_Image_Double(I(n+1:256*256)) = [];

N = 2*n+2000;
m = 2000;
L = [1:round(percenterasures*N)];
LC = setdiff([1:N],L);
Times = zeros(1,5);
Errors = zeros(1,7);

f = Compressed_Image_Double;

fprintf('Creating Frames \n');

A = randn(N,n+m);
[A,~] = qr(A,0);

F = sqrt(N/n)*A(:,1:n)';
G = (n/N)*DF;
M = sqrt(N/m)*A(:,n+1:n+m)';

fprintf('Reconstructing Erasures \n');

FC = G' * f;
FC(L) = zeros(size(L'));
noiselessf_R = F * FC;
Errors(1) = norm(f-noiselessf_R);
noise = randn(size(LC'));
noise = snr * noise ./ norm(noise) * norm(FC(LC));
FC(LC) = FC(LC) + noise;
FC_NDB = FC;
FC_RDI = FC;
FC_RDIN = FC;
FC_ERM = FC;
FC_FORC = FC;
f_R = F*FC;
Errors(2) = norm(f-f_R);

% Nilpotent Double Bridging Reconstruction

tic;
FRCL = G(:,L)' * f_R;
FRCB = G(:,W)' * f_R;
C_NDB = pinv(F(:,L)'*G(:,W))*(F(:,L)'*G(:,L));
FC_NDB(L) = C_NDB' * (FC_NDB(W) - FRCB) + FRCL;
g_NDB = f_R + F(:,L) * FC_NDB(L);
Times(1) = toc;
Errors(3) = norm(f-g_NDB);

% Reduced Direct Inversion Reconstruction

tic;
M_RDI = G(:,L)' * F(:,L);
C_RDI = (eye(length(L)) - M_RDI) \ eye(length(L));
g_RDI = f_R + F(:,L) * (C_RDI * (G(:,L)' * f_R));
Times(2) = toc;
Errors(4) = norm(f-g_RDI);

% Reduced Direct Inversion with Neumann Iterations

tic;
M_RDIN = G(:,L)' * F(:,L);
Mnorm = norm(M_RDIN);
NumIter = round(log(tolerance*(1-Mnorm))/log(Mnorm));
C_m = eye(length(L));
for(j = 1:1:NumIter)
    C_m = eye(length(L)) + M_RDIN * C_m;
end
g_RDIN = f_R + F(:,L) * (C_m * (G(:,L)' * f_R));
Times(3) = toc;
Errors(5) = norm(f-g_RDIN)


% Erasure Recovery Matrices Reconstruction

tic;
FC_ERM(L) = -(M(:,L)' * M(:,L))\(M(:,L)' * (M(:,LC) * FC(LC)));
g_ERM = f_R + F(:,L) * FC(L);
Times(4) = toc;
Errors(6) = norm(f-g_ERM);

% FORC Method Reconstruction

tic;
g_FORC = (G(:,LC) * G(:,LC)') \ (G(:,LC) * FC_FORC)
Times(5) = toc;
Errors(7) = norm(f-g_FORC);

fprintf('Plotting Images \n');

I1 = sort(I(1:n),'ascend');

C_f = zeros(256*256,1); % Compressed Image.
C_f(I1) = f;
Uncompressed_f = ifft(C_f);
Uncompressed_f = reshape(Uncompressed_f,[256,256]);
J_f = uint8(Uncompressed_f);

C_noiselessf_R = zeros(256*256,1);
C_noiselessf_R(I1(1:n)) = noiselessf_R;
Uncompressed_noiselessf_R = ifft(C_noiselessf_R);
Uncompressed_noiselessf_R = reshape(Uncompressed_noiselessf_R,[256,256]);
J_noiselessf_R = uint8(Uncompressed_noiselessf_R);

C_f_R = zeros(256*256,1); % Compressed Image.
C_f_R(I1) = f_R;
Uncompressed_f_R = ifft(C_f_R);
Uncompressed_f_R = reshape(Uncompressed_f_R,[256,256]);
J_f_R = uint8(Uncompressed_f_R);

C_g_NDB = zeros(256*256,1); % Compressed Image.
C_g_NDB(I1) = g_NDB;
Uncompressed_g_NDB = ifft(C_g_NDB);
Uncompressed_g_NDB = reshape(Uncompressed_g_NDB,[256,256]);
J_g_NDB = uint8(Uncompressed_g_NDB);

C_g_RDI = zeros(256*256,1); % Compressed Image.
C_g_RDI(I1) = g_RDI;
Uncompressed_g_RDI = ifft(C_g_RDI);
Uncompressed_g_RDI = reshape(Uncompressed_g_RDI,[256,256]);
J_g_RDI = uint8(Uncompressed_g_RDI);

C_g_RDIN = zeros(256*256,1); % Compressed Image.
C_g_RDIN(I1) = g_RDIN;
Uncompressed_g_RDIN = ifft(C_g_RDIN);
Uncompressed_g_RDIN = reshape(Uncompressed_g_RDIN,[256,256]);
J_g_RDIN = uint8(Uncompressed_g_RDIN);

C_g_ERM = zeros(256*256,1); % Compressed Image.
C_g_ERM(I1) = g_ERM;
Uncompressed_g_ERM = ifft(C_g_ERM);
Uncompressed_g_ERM = reshape(Uncompressed_g_ERM,[256,256]);
J_g_ERM = uint8(Uncompressed_g_ERM);

C_g_FORC = zeros(256*256,1); % Compressed Image.
C_g_FORC(I1) = g_FORC;
Uncompressed_g_FORC = ifft(C_g_FORC);
Uncompressed_g_FORC = reshape(Uncompressed_g_FORC,[256,256]);
J_g_FORC = uint8(Uncompressed_g_FORC);

figure;

subplot(2,5,2);
imshow(J_f);
title('Compressed Image');

subplot(2,5,3);
imshow(J_noiselessf_R);
title('Noise Free Partial Reconstruction');

subplot(2,5,4);
imshow(J_f_R);
title('Noisy Partial Reconstruction');

subplot(2,5,6);
imshow(J_g_NDB);
title('Nilpotent Double Bridging');

subplot(2,5,7);
imshow(J_g_RDI);
title('Reduced Direct Inversion');

subplot(2,5,8);
imshow(J_g_RDIN);
title('Neumann Reduced Direct Inversion');

subplot(2,5,9);
imshow(J_g_ERM);
title('Erasure Recovery Matrices');

subplot(2,5,10);
imshow(J_g_FORC);
title('FORC Method');






