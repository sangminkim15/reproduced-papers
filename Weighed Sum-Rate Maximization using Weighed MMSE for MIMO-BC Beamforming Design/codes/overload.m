function [] = overload ()

% Overloaded Case %

SNRdB = -10 : 5 : 30;
SNR = 10.^(SNRdB / 10);
Etx = 1;
WSR = zeros(1000, length(SNR));
WSRTxMF = zeros(1000, length(SNR));

u = ones(1, 20);

for i = 1 : length(SNR)
    for j = 1 : 1000
        H = channel(4, 1, 20, SNR(i));
        WSR(j,i) = WSRBF_WMMSE2 (u, H, Etx);
        WSRTxMF(j,i) = TxMF (u, H, Etx);
    end
end

WSRmean = mean(WSR);
WSRTxMean = mean(WSRTxMF);

plot(SNRdB, WSRmean, '-*', SNRdB, WSRTxMean, '--o', 'LineWidth', 2);
title('Overloaded Case');
xlabel('SNR (dB)');
ylabel('Sum-Rate (bits/complex dim)');
legend('WSRBF-WMMSE2 (10 iterations - TxMF init)', 'Transmit Matched Filter (TxMF)');
ylim([0 45]);
grid on;

end