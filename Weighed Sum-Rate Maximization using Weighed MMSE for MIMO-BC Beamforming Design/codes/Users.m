function [] = Users ()

% Performance w.r.t. Number of Users %

SNRdB = -10 : 5 : 30;
SNR = 10.^(SNRdB / 10);
Etx = 1;

WSR11 = zeros(1000, length(SNR));
WSR12 = zeros(1000, length(SNR));
WSR16 = zeros(1000, length(SNR));

WSR21 = zeros(1000, length(SNR));
WSR22 = zeros(1000, length(SNR));
WSR26 = zeros(1000, length(SNR));

u1 = ones(1, 1);
u2 = ones(1, 2);
u6 = ones(1, 6);

for i = 1 : length(SNR)
    for j = 1 : 1000
        H1 = channel(4, 2, 1, SNR(i));
        H2 = channel(4, 2, 2, SNR(i));
        H6 = channel(4, 2, 6, SNR(i));
        
        WSR11(j,i) = WSRBF_WMMSE1 (u1, H1, Etx);
        WSR12(j,i) = WSRBF_WMMSE1 (u2, H2, Etx);
        WSR16(j,i) = WSRBF_WMMSE1 (u6, H6, Etx);
        
        WSR21(j,i) = WSRBF_WMMSE2 (u1, H1, Etx);
        WSR22(j,i) = WSRBF_WMMSE2 (u2, H2, Etx);
        WSR26(j,i) = WSRBF_WMMSE2 (u6, H6, Etx);
    end
end

WSRmean11 = mean(WSR11);
WSRmean12 = mean(WSR12);
WSRmean16 = mean(WSR16);

WSRmean21 = mean(WSR21);
WSRmean22 = mean(WSR22);
WSRmean26 = mean(WSR26);

plot(SNRdB, WSRmean11, 'r-*', SNRdB, WSRmean12, 'r:', SNRdB, WSRmean16, 'r--o', SNRdB, WSRmean21, 'b-*', SNRdB, WSRmean22, 'b:', SNRdB, WSRmean26, 'b--o', 'LineWidth', 2);
title('Performance w.r.t. Number of Users');
xlabel('SNR (dB)');
ylabel('Sum-Rate (bits/complex dim)');
legend('WSRBF-WMMSE1, K=1', 'WSRBF-WMMSE1, K=2', 'WSRBF-WMMSE1, K=6', 'WSRBF-WMMSE2, K=1', 'WSRBF-WMMSE2, K=2', 'WSRBF-WMMSE2, K=6');
ylim([0 45]);
grid on;

end