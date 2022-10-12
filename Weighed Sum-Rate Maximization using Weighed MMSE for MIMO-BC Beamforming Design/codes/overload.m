function [] = overload ()

% Overloaded Case %

warning off

SNRdB = -10 : 5 : 30;
SNR = 10.^(SNRdB / 10);
Etx = 1;
WSR = zeros(1000, length(SNR));
WSRTxMF = zeros(1000, length(SNR));
WSRTxZF = zeros(1000, length(SNR));

u = ones(1, 10);

for i = 1 : length(SNR)
    for j = 1 : 1000
        H = channel(8, 1, 10, SNR(i));
        WSR(j,i) = WSRBF_WMMSE2 (u, H, Etx);
        WSRTxMF(j,i) = TxMF (u, H, Etx);
        WSRTxZF(j,i) = TxZF (u, H, Etx);
    end
end

WSRmean = mean(WSR);
WSRMFMean = mean(WSRTxMF);
WSRZFMean = mean(WSRTxZF);

plot(SNRdB, WSRmean, '-*', SNRdB, WSRMFMean, '--o', SNRdB, WSRZFMean, 'LineWidth', 2);
title('Overloaded Case');
xlabel('SNR (dB)');
ylabel('Sum-Rate (bits/complex dim)');
legend('WSRBF-WMMSE2 (10 iterations - TxMF init)', 'Transmit Matched Filter (TxMF)', 'Transmit Zero-Forcing Filter (TxZF)');
grid on;

end