function [] = convergence ()

iterations = cell(100, 1);
Zplot = cell(100, 1);
SNRplot = cell(100, 1);
SNRmaxplot = cell(100, 1);
SNReqplot = cell(100, 1);
ERRORplot = cell(100, 1);
ERRORplot2 = cell(100, 1);
SS = cell(100, 1);
iterations_len = zeros(100, 1);

for idx = 1 : 100
    [iterations{idx}, Zplot{idx}, SNRplot{idx}, SNRmaxplot{idx}, SNReqplot{idx}, ERRORplot{idx}, ERRORplot2{idx}, SS{idx}] = cyclic ();
    iterations_len(idx, 1) = length(iterations{idx});
end

min_len = min(iterations_len);
it = zeros(1, min_len);

for jdx = 1 : min_len
    it(jdx) = jdx;
end

Zpl = zeros(100, min_len);
SNRpl = zeros(100, min_len);
SNRmaxpl = zeros(100, min_len);
SNReqpl = zeros(100, min_len);
ERRORpl = zeros(100, min_len);
ERRORpl2 = zeros(100, min_len);

for kdx = 1 : 100
    ZZpl = Zplot{kdx, 1};
    SSNRpl = SNRplot{kdx, 1};
    SSNRmaxpl = SNRmaxplot{kdx, 1};
    SSNReqpl = SNReqplot{kdx, 1};
    EERRORpl = ERRORplot{kdx, 1};
    EERRORpl2 = ERRORplot2{kdx, 1};
    
    Zpl(kdx, :) = ZZpl(1, 1:min_len);
    SNRpl(kdx, :) = SSNRpl(1, 1:min_len);
    SNRmaxpl(kdx, :) = SSNRmaxpl(1, 1:min_len);
    SNReqpl(kdx, :) = SSNReqpl(1, 1:min_len);
    ERRORpl(kdx, :) = EERRORpl(1, 1:min_len);
    ERRORpl2(kdx, :) = EERRORpl2(1, 1:min_len);

end

Zmean = mean(Zpl);
SNRmean = mean(SNRpl);
SNRmaxmean = mean(SNRmaxpl);
SNReqmean = mean(SNReqpl);
ERRORmean = mean(ERRORpl);
ERRORmean2 = mean(ERRORpl2);

figure
plot(it, Zmean, 'LineWidth', 1.5);
xlabel('# of iterations');
ylabel('Objective Function');
title('Objective Function Convergence');
grid on

figure
plot(it, SNRmean, it, SNRmaxmean, it, SNReqmean, 'LineWidth', 1.5);
xlabel('# of iterations');
ylabel('Radar SNR (dB)');
title('Radar SNR Convergence');
legend('Proposed Algorithm', 'Radar SNR Upper Bound', 'Equal Power per Subcarrier', 'location', 'southeast');
xlim([0 50]);
ylim([28 32]);
grid on

figure
semilogy(it, ERRORmean, it, ERRORmean2, 'LineWidth', 1.5);
xlabel('# of iterations');
ylabel('Error Probability');
title('Error Probability Convergence');
legend('Proposed Algorithm', 'Equal Power per Subcarrier');
xlim([0 200]);
ylim([1e-4 2e-1]);
grid on

figure
[xx, yy] = meshgrid(1:64, 1:64);
plot3(xx, yy, abs(SS{1}), 'LineWidth', 1.5);
zlabel('|{S^H}S|')
title('Diagonality of {S^H}S');
grid on
    
end