function [] = threshold ()

SNRdB = [-30 -20 -10 0];
SNRmin = 10 .^ (SNRdB ./ 10);

SNRrad = zeros(100, length(SNRmin), 500);
ERRORcomm = zeros(100, length(SNRmin), 500);

for idx = 1 : length(SNRmin)
    for jdx = 1 : 100
        [~, ~, SNR_rad, ERROR_comm, ~] = cyclic3 (SNRmin(idx));
        for kdx = 1 : 500
            SNRrad(jdx, idx, kdx) = SNR_rad(kdx);
            ERRORcomm(jdx, idx, kdx) = ERROR_comm(kdx);
        end
    end
end

SNRradplot = mean(SNRrad);
ERRORcommplot = mean(ERRORcomm);

iterations = zeros(1, 500);
for ldx = 1 : 500
    iterations(ldx) = ldx;
end

SNRrad1 = zeros(1, 500);
SNRrad2 = zeros(1, 500);
SNRrad3 = zeros(1, 500);
SNRrad4 = zeros(1, 500);
 
ERRORcomm1 = zeros(1, 500);
ERRORcomm2 = zeros(1, 500);
ERRORcomm3 = zeros(1, 500);
ERRORcomm4 = zeros(1, 500);

for iidx = 1 : 500
    SNRrad1(iidx) = SNRradplot(1, 1, iidx);
    SNRrad2(iidx) = SNRradplot(1, 2, iidx);
    SNRrad3(iidx) = SNRradplot(1, 3, iidx);
    SNRrad4(iidx) = SNRradplot(1, 4, iidx);
    
    
    ERRORcomm1(iidx) = ERRORcommplot(1, 1, iidx);
    ERRORcomm2(iidx) = ERRORcommplot(1, 2, iidx);
    ERRORcomm3(iidx) = ERRORcommplot(1, 3, iidx);
    ERRORcomm4(iidx) = ERRORcommplot(1, 4, iidx);
end

figure
plot(iterations, SNRrad1, iterations, SNRrad2, iterations, SNRrad3, iterations, SNRrad4, 'LineWidth', 1.5);
xlabel('iterations');
ylabel('Radar SNR (dB)');
legend('SNR_{min} = -30dB', 'SNR_{min} = -20dB', 'SNR_{min} = -10dB', 'SNR_{min} = -0dB');
title('Radar SNR Performance');
grid on;

figure
semilogy(iterations, ERRORcomm1, iterations, ERRORcomm2, iterations, ERRORcomm3, iterations, ERRORcomm4, 'LineWidth', 1.5);
xlabel('iterations)');
ylabel('BER');
legend('SNR_{min} = -30dB', 'SNR_{min} = -20dB', 'SNR_{min} = -10dB', 'SNR_{min} = -0dB');
title('Communication BER Performance');
grid on;

end