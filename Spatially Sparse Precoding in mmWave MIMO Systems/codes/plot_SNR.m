function [] = plot_SNR ()

warning off;

% SNR
rhodB = -40 : 4 : 0;

[Iss1, Iun1] = SNR (1, 64, 16, 4, 4, 8, 10, rhodB);         % Ns = 1, 64 * 16, 4 RF chains
[Iss2, Iun2] = SNR (2, 64, 16, 4, 4, 8, 10, rhodB);         % Ns = 2, 64 * 16, 4 RF chains

Xss1 = plot(rhodB, Iss1, 'r');
Xss1.LineWidth = 2;
Xss1.Marker = 'o';
Xss1.MarkerEdgeColor = 'r';
hold on;

Xun1 = plot(rhodB, Iun1, 'b');
Xun1.LineWidth = 2;
Xun1.Marker = 's';
Xun1.MarkerEdgeColor = 'b';

Xss2 = plot(rhodB, Iss2, 'r');
Xss2.LineWidth = 2;
Xss2.Marker = 'o';
Xss2.MarkerEdgeColor = 'r';

Xun2 = plot(rhodB, Iun2, 'b');
Xun2.LineWidth = 2;
Xun2.Marker = 's';
Xun2.MarkerEdgeColor = 'b';
hold off;

legend([Xss1 Xun1], 'SS Precoding & Combining', 'Optimal Unstrained Precoding', 'Location', 'northwest');
title('Spatially Sparse Precoding  vs Unconstrained Optimum Precoding');
xlabel('SNR (dB)');
ylabel('Spectral Efficiency (bits/s/Hz)');
grid on;

end