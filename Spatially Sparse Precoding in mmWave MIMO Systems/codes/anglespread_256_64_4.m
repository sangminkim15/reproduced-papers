function [] = anglespread_256_64_4()

warning off;

% # of antennas
Nt = 256;
Nr = 64;
sqrtNt = sqrt(Nt);
sqrtNr = sqrt(Nr);

% # of RF chains
NtRF = 4;
NrRF = 4;

% # of departed(arrived) rays
Ncl = 8;
Nray = 10;

% Angle Spread (std.)
std = pi/180 : pi/180 : pi/12;
stddeg = std.* (180/pi);

% Noise
sigma_n = 1;    % 0dB

rho = 1;

I1 = zeros(1000, length(std));     % SS (Spatially Sparse) / Ns = 1
I2 = zeros(1000, length(std));     % SS (Spatially Sparse) / Ns = 2
I3 = zeros(1000, length(std));     % Unconstrained / Ns = 1
I4 = zeros(1000, length(std));     % Unconstrained / Ns = 2

for i = 1 : length(std)
    for l = 1 : 1000
        Atcell =ArrayResponse_cell(sqrtNt, Ncl, Nray, std(i));
        Arcell =ArrayResponse_cell(sqrtNr, Ncl, Nray, std(i));
        
        At = cell2mat(Atcell);
        Ar = cell2mat(Arcell);
        
        % CHANNEL Formation (Lines 39 ~ 51)
        H = zeros(Nr, Nt);
        
        for p = 1 : Ncl
            Atmat = Atcell{1,p};
            Armat = Arcell{1,p};
    
            for q = 1 : Nray
                alpha = sqrt(1/2) * (randn(1,1) + 1i * randn(1,1));
        
                H = H + alpha * Armat(:,q) * Atmat(:,q)';
            end
        
        end
        
        H = (sqrt(Nt * Nr) / norm(H, 'fro')) * H;
        
        % SS - Spatial Sparse Precoding / Decoding
        I1(l,i) = SS(1, NtRF, NrRF, H, At, Ar, sigma_n, rho);
        I2(l,i) = SS(2, NtRF, NrRF, H, At, Ar, sigma_n, rho);
        
        % unconstrained - Unconstrained Precoding / Decoding
        I3(l,i) = unconstrained (1, H, sigma_n, rho);
        I4(l,i) = unconstrained (2, H, sigma_n, rho);    
    end
end

I1mean = mean(I1);
I2mean = mean(I2);
I3mean = mean(I3);
I4mean = mean(I4);

x = plot(stddeg, I1mean, stddeg, I2mean, stddeg, I3mean, stddeg, I4mean, 'b');

x(1).LineWidth = 2;
x(2).LineWidth = 2;
x(3).LineWidth = 2;
x(4).LineWidth = 2;

x(1).Marker = 'o';
x(2).Marker = 'o';
x(3).Marker = 's';
x(4).Marker = 's';

legend('SS Precoding & Combining, Ns=1', 'SS Precoding & Combining, Ns=2', 'Optimal Unstrained Precoding, Ns=1', 'Optimal Unstrained Precoding, Ns=2');
title('Spatially Sparse Precoding  vs Unconstrained Optimum Precoding');
xlabel('Angle Spread (\circ)');
ylabel('Spectral Efficiency (bits/s/Hz)');
grid on;

end