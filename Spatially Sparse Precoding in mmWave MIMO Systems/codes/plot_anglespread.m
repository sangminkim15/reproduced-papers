function [] = plot_anglespread ()

warning off;

% Angle Spread (std.)
std = pi/180 : pi/180 : pi/12;
stddeg = std.* (180/pi);

[Iss1, Iun1] = anglespread (1, 64, 16, 4, 4, 8, 10, std);       % Ns = 1, 64 * 16, 4 RF chains
[Iss2, Iun2] = anglespread (2, 64, 16, 4, 4, 8, 10, std);       % Ns = 2, 64 * 16, 4 RF chains
[Iss3, Iun3] = anglespread (1, 256, 64, 4, 4, 8, 10, std);      % Ns = 1, 256 * 64, 4 RF chains
[Iss4, ~]    = anglespread (2, 64, 16, 6, 6, 8, 10, std);       % Ns = 2, 64 * 16, 6 RF chains

Xss1 = plot(stddeg, Iss1, 'r');
Xss1.LineWidth = 2;
Xss1.Marker = 'o';
Xss1.MarkerEdgeColor = 'r';
hold on;

Xun1 = plot(stddeg, Iun1, 'b');
Xun1.LineWidth = 2;
Xun1.Marker = 's';
Xun1.MarkerEdgeColor = 'b';

Xss2 = plot(stddeg, Iss2, 'r');
Xss2.LineWidth = 2;
Xss2.Marker = 'o';
Xss2.MarkerEdgeColor = 'r';

Xun2 = plot(stddeg, Iun2, 'b');
Xun2.LineWidth = 2;
Xun2.Marker = 's';
Xun2.MarkerEdgeColor = 'b';

Xss3 = plot(stddeg, Iss3, 'r');
Xss3.LineWidth = 2;
Xss3.Marker = 'o';
Xss3.MarkerEdgeColor = 'r';

Xun3 = plot(stddeg, Iun3, 'b');
Xun3.LineWidth = 2;
Xun3.Marker = 's';
Xun3.MarkerEdgeColor = 'b';

Xss4 = plot(stddeg, Iss4, 'r');
Xss4.LineWidth = 2;
Xss4.Marker = 'o';
Xss4.MarkerEdgeColor = 'r';
hold off;

legend([Xss1 Xun1], 'SS Precoding & Combining', 'Optimal Unstrained Precoding', 'Location', 'southwest');
title('Spatially Sparse Precoding  vs Unconstrained Optimum Precoding');
xlabel('Angle Spread (\circ)');
ylabel('Spectral Efficiency (bits/s/Hz)');
grid on;

end