function [] = fullload1 ()

% Fully-Loaded Case %
% For Highly-Correlated Channels (UE1 - UE4) %

SNRdB = -10 : 5 : 30;
SNR = 10.^(SNRdB / 10);
Etx = 1;
WSR1 = zeros(1000, length(SNR));
WSR2 = zeros(1000, length(SNR));
WSRTxMF = zeros(1000, length(SNR));
WSRTxZF = zeros(1000, length(SNR));

u = ones(1, 4);

for idx = 1 : length(SNR)
    for jdx = 1 : 1000
        H13 = channel(64, 1, 3, SNR(idx));
        var = norm(H13{1,1}) * 0.01;
        H4 = H13{1,1} + sqrt(var/2)*(randn(1, 64)+1i*randn(1, 64));
        
        H = cell(4, 1);
        for iidx = 1 : 3
            H{iidx,1} = H13{iidx, 1};
        end
        H{4,1} = H4;
        
        WSR1(jdx,idx) = WSRBF_WMMSE1 (u, H, Etx);
        WSR2(jdx,idx) = WSRBF_WMMSE2 (u, H, Etx);
        WSRTxMF(jdx,idx) = TxMF (u, H, Etx);
        WSRTxZF(jdx,idx) = TxZF (u, H, Etx);
    end
end

WSRmean1 = mean(WSR1);
WSRmean2 = mean(WSR2);
WSRMFMean = mean(WSRTxMF);
WSRZFMean = mean(WSRTxZF);

plot(SNRdB, WSRmean1,'-*', SNRdB, WSRmean2, '-*', SNRdB, WSRMFMean, '--o', SNRdB, WSRZFMean, '--v', 'LineWidth', 2);
title('Fully-Loaded Case');
xlabel('SNR (dB)');
ylabel('Sum-Rate (bits/complex dim)');
legend('WSRBF-WMMSE1 (convergence / 10 random init)', 'WSRBF-WMMSE2 (10 iterations - TxMF init)', 'Transmit Matched Filter (TxMF)', 'Transmit Zero-Forcing Filter (TxZF)');
grid on;

end