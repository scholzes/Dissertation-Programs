m = 250;
n = 250;
N = 1000;
snr = .05;
Trials = 50;
EC = [10:10:240];

TimeDataERM = zeros(Trials,length(EC));
TimeDataB2 = zeros(Trials,length(EC));
TimeDataRDI = zeros(Trials,length(EC));
TimeDataRDIN = zeros(Trials,length(EC));
TimeDataFORC = zeros(Trials,length(EC));
DataERM = zeros(Trials,length(EC));
DataB2 = zeros(Trials,length(EC));
DataRDI = zeros(Trials,length(EC));
DataRDIN = zeros(Trials,length(EC));
DataFORC = zeros(Trials,length(EC));

for(k=1:1:length(EC))
    
    for(t = 1:1:Trials)
        
        L = [1:1:EC(k)];
        LC = setdiff(1:N,L);
        
        % ERM
        
        M = (1/sqrt(m)) * randn(m,N);
        A = [M',randn(N,n)];
        [A,~] = qr(A,0);
        F1 = A(:,m+1:m+n)';

        f = randn(n,1);
        f = f./norm(f);
        
        FC1 = F1' * f;
        noise = randn(length(LC),1);
        noise = snr * norm(FC1(LC)) / norm(noise) * noise;
        FC1(LC) = FC1(LC) + noise;
        FC1(L) = zeros(size(L'));
        f_R = F1 * FC1;
        
        tic;
        LC = setdiff(1:N,L);
        FC1(L) = -(M(:,L)' * M(:,L))\(M(:,L)' * (M(:,LC) * FC1(LC)));
        gERM = f_R + F1(:,L) * FC1(L);
        TimeDataERM(t,k) = toc;
        DataERM(t,k) = norm(f-gERM);
        
        % Other Methods
        
        F = randn(N,n);
        [F,~] = qr(F,0);
        F = sqrt(N/n)*F(:,1:n)';
        G = (n/N)*F;
        
        FC = G' * f;
        noise = randn(size(LC'));
        noise = snr * norm(FC(LC))/norm(noise) * noise;
        FC(LC) = FC(LC) + noise;
        FC(L) = zeros(size(L'));
        FCB2 = FC;
        FCRDI = FC;
        FCRDIN = FC;
        FCFORC = FC;
        f_R = F*FC;
        
        % 25% Overbridging
        tic;
        W = [EC(k)+1:1:3*EC(k)];
        FRCL = G(:,L)' * f_R;
        FRCB = G(:,W)' * f_R;
        CB2 = (F(:,L)'*G(:,W))\(F(:,L)'*G(:,L));
        FCB2(L) = CB2' * (FC(W) - FRCB) + FRCL;
        gB2 = f_R + F(:,L) * FCB2(L);
        TimeDataB2(t,k) = toc;
        DataB2(t,k) = norm(f-gB2);
        
        % RDI
        tic;
        MRDI = G(:,L)' * F(:,L);
        gRDI = f_R + F(:,L) * ((eye(length(L)) - MRDI) \ (G(:,L)' * f_R));
        TimeDataRDI(t,k) = toc;
        DataRDI(t,k) = norm(f-gRDI);
        
        % RDIN
        tic;
        tolerance = 10^(-5);
        MRDIN = G(:,L)' * F(:,L);
        Mnorm = max(abs(eigs(MRDIN)));
        NumIter = round(log(tolerance*(1-Mnorm))/log(Mnorm));
        g0 = G(:,L)' * f_R;
        Cg_RDIN = zeros(size(L'));
        for(j = 1:1:NumIter)
            Cg_RDIN = g0 + MRDIN * Cg_RDIN;
        end
        gRDIN = f_R + F(:,L) * Cg_RDIN;
        TimeDataRDIN(t,k) = toc;
        DataRDIN(t,k) = norm(f-gRDIN);
        
        % FORC
        tic;
        gFORC = (G(:,LC) * G(:,LC)') \ (G(:,LC) * FCFORC(LC));
        TimeDataFORC(t,k) = toc;
        DataFORC(t,k) = norm(f-gFORC);

    end
    
    k

end

figure;
plot(EC,mean(DataERM),'*');
hold on;
plot(EC,mean(DataB2),'^');
plot(EC,mean(DataRDI),'o');
plot(EC,mean(DataRDIN),'+');
plot(EC,mean(DataFORC),'x');
title('Erasure Set Size vs Reconstruction Error');
xlabel('Erasure Set Size');
ylabel('Reconstruction Error');
legend('Erasure Recovery Matrices','25 Percent Overbridging','Reduced Direct Inversion','Reduced Direct Inversion with Neumann','FORC Method','Location','northwest')
hold off;

figure;
plot(EC,mean(TimeDataERM),'*');
hold on;
plot(EC,mean(TimeDataB2),'^');
plot(EC,mean(TimeDataRDI),'o');
plot(EC,mean(TimeDataRDIN),'+');
plot(EC,mean(TimeDataFORC),'x');
title('Erasure Set Size vs Reconstruction Time');
xlabel('Erasure Set Size');
ylabel('Time (in Seconds)');
legend('Erasure Recovery Matrices','25 Percent Overbridging','Reduced Direct Inversion','Reduced Direct Inversion with Neumann','FORC Method','Location','northwest')
hold off;

