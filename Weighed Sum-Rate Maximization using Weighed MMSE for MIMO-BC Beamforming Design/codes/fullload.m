function [] = fullload ()

% Fully-Loaded Case %

SNRdB = -10 : 5 : 30;
SNR = 10.^(SNRdB / 10);
Etx = 1;
WSR1 = zeros(1000, length(SNR));
WSR2 = zeros(1000, length(SNR));
WSRTxMF = zeros(1000, length(SNR));

u = ones(1, 4);

for i = 1 : length(SNR)
    for j = 1 : 1000
        H = channel(4, 1, 4, SNR(i));
        WSR1(j,i) = WSRBF_WMMSE1 (u, H, Etx);
        WSR2(j,i) = WSRBF_WMMSE2 (u, H, Etx);
        WSRTxMF(j,i) = TxMF (u, H, Etx);
    end
end

WSRmean1 = mean(WSR1);
WSRmean2 = mean(WSR2);
WSRTxMean = mean(WSRTxMF);

plot(SNRdB, WSRmean1,'-*', SNRdB, WSRmean2, '-*', SNRdB, WSRTxMean, '--o', 'LineWidth', 2);
title('Fully-Loaded Case');
xlabel('SNR (dB)');
ylabel('Sum-Rate (bits/complex dim)');
legend('WSRBF-WMMSE1 (convergence / 10 random init)', 'WSRBF-WMMSE2 (10 iterations - TxMF init)', 'Transmit Matched Filter (TxMF)');
ylim([0 40]);
grid on;

end