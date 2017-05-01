clc;

% Original Images are 256 pixels X 256 pixels.

fprintf('Reading Image \n');

COMPRESSION_PERCENT = 0.15; % Compressed Signal will be approximately
% n = 256^2 * COMPRESSION_PERCENT dimensional.
snr = .05;
percenterasures = .05;

Original_Image_Double = double(imread('Pepper.bmp'));

fprintf('Performing Image Compression \n');

Compressed_Image_Double = fft(reshape(Original_Image_Double,[256*256,1]));
[S,I] = sort(abs(Compressed_Image_Double),'descend');
n = round(COMPRESSION_PERCENT*256*256);
LSC = Compressed_Image_Double(I(n+1:256*256));
Compressed_Image_Double(I(n+1:256*256)) = [];

m = 2000;
N = 2*n+m;
L = [1:round(percenterasures*N)];
LC = setdiff([1:N],L);
Times = zeros(1,5);
Errors = zeros(1,8);

f = Compressed_Image_Double;

fprintf('Creating Frames \n');

A = randn(N,n+m);
[A,~] = qr(A,0);

F = sqrt(N/n)*A(:,1:n)';
G = (n/N)*F;
M = sqrt(N/m)*A(:,n+1:n+m)';

fprintf('Reconstructing Erasures \n');

FC = G' * f;
FCNonErased = FC;
FC(L) = zeros(size(L'));
noiselessf_R = F * FC;
Errors(1) = norm(f-noiselessf_R);
noise = randn(size(LC'));
noise = snr * noise ./ norm(noise) * norm(FC(LC));
FCNonErased(LC) = FCNonErased(LC) + noise;
g_NonErased = F * FCNonErased;
Errors(8) = norm(f - g_NonErased);
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
W = [length(L)+1:3*length(L)];
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
% C_RDI = (eye(length(L)) - M_RDI) \ eye(length(L));
g_RDI = f_R + F(:,L) * ((eye(length(L)) - M_RDI) \ (G(:,L)' * f_R));
Times(2) = toc;
Errors(4) = norm(f-g_RDI);

% Reduced Direct Inversion with Neumann Iterations

tic;
tolerance = 10^(-5);
M_RDIN = G(:,L)' * F(:,L);
toc;
Mnorm = max(abs(eigs(M_RDIN)));
toc;
NumIter = round(log(tolerance*(1-Mnorm))/log(Mnorm));
g0 = G(:,L)' * f_R;
Cg_RDIN = zeros(size(L'));
for(j = 1:1:NumIter)
    Cg_RDIN = g0 + M_RDIN * Cg_RDIN;
end
g_RDIN = f_R + F(:,L) * Cg_RDIN;
Times(3) = toc;
Errors(5) = norm(f-g_RDIN);


% Erasure Recovery Matrices Reconstruction

tic;
FC_ERM(L) = -(M(:,L)' * M(:,L))\(M(:,L)' * (M(:,LC) * FC_ERM(LC)));
g_ERM = f_R + F(:,L) * FC_ERM(L);
Times(4) = toc;
Errors(6) = norm(f-g_ERM);

% FORC Method Reconstruction

tic;
g_FORC = (G(:,LC) * G(:,LC)') \ (G(:,LC) * FC_FORC(LC));
Times(5) = toc;
Errors(7) = norm(f-g_FORC);

fprintf('Plotting Images \n');

I1 = sort(I(1:n),'ascend');

C_f = zeros(256*256,1); % Compressed Image.
C_f(I1) = f;
Uncompressed_f = ifft(C_f);
Uncompressed_f = reshape(Uncompressed_f,[256,256]);
J_f = uint8(Uncompressed_f);

C_g_NonErased = zeros(256*256,1); % Compressed Image.
C_g_NonErased(I1) = g_NonErased;
Uncompressed_g_NonErased = ifft(C_g_NonErased);
Uncompressed_g_NonErased = reshape(Uncompressed_g_NonErased,[256,256]);
J_g_NonErased = uint8(Uncompressed_g_NonErased);

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

subplot(2,5,1);
imshow(J_f);
title('Compressed Image');

subplot(2,5,3);
imshow(J_g_NonErased);
title('Noise Only');

subplot(2,5,4);
imshow(J_noiselessf_R);
title('Noise Free Partial Reconstruction');

subplot(2,5,5);
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

Errors = Errors ./ norm(f)
Times