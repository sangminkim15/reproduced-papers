function [] = performance ()

SNRdB = 0 : 2 : 10;
SNRin = 10 .^ (SNRdB ./ 10);

SNRraddB = zeros(100, length(SNRin));
SNRmaxdB = zeros(100, length(SNRin));
SNReqdB = zeros(100, length(SNRin));
ErrorP = zeros(100, length(SNRin));
ErrorP2 = zeros(100, length(SNRin));

for idx = 1 : length(SNRin)
    for jdx = 1 : 100
        [~, SNRraddB(jdx, idx), SNRmaxdB(jdx, idx), SNReqdB(jdx, idx), ~, ErrorP(jdx, idx), ErrorP2(jdx, idx)] = cyclic2 (SNRin(idx));
    end
end

SNRradplot = mean(SNRraddB);
SNRmaxplot = mean(SNRmaxdB);
SNReqplot = mean(SNReqdB);
ErrorPplot = mean(ErrorP);
ErrorP2plot = mean(ErrorP2);

figure
plot(SNRdB, SNRradplot, '-x', SNRdB, SNRmaxplot, '-o', SNRdB, SNReqplot, 'LineWidth', 1.5);
xlabel('input SNR (dB)');
ylabel('Radar SNR (dB)');
legend('Proposed Algorithm', 'Radar SNR Upper Bound', 'Equal Power per Subcarrier', 'Location', 'northwest');
title('Radar SNR Performance');
grid on;

figure
semilogy(SNRdB, ErrorPplot, '-x', SNRdB, ErrorP2plot, '-o', 'LineWidth', 1.5);
xlabel('input SNR (dB)');
ylabel('BER');
legend('Proposed Algorithm', 'Equal Power per Subcarrier', 'Location', 'southwest');
title('Communication BER Performance');
grid on;

end