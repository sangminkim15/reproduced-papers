function [] = SNR_64_16_4()

% # of antennas
Nt = 64;
Nr = 16;
sqrtNt = sqrt(Nt);
sqrtNr = sqrt(Nr);

% # of RF chains
NtRF = 4;
NrRF = 4;

% # of departed(arrived) rays
Ncl = 8;
Nray = 10;

% Angle Spread (std.)
std = pi/24;

% SNR
sigma_n = 1;                       % 0dB
rhodB = -40 : 4 : 0;
rho = 10.^(rhodB/10);

I1 = zeros(1000, length(rho));     % SS (Spatially Sparse) / Ns = 1
I2 = zeros(1000, length(rho));     % SS (Spatially Sparse) / Ns = 2
I3 = zeros(1000, length(rho));     % Unconstrained / Ns = 1
I4 = zeros(1000, length(rho));     % Unconstrained / Ns = 2

for i = 1 : length(rho)
    for l = 1 : 1000
        Atcell = ArrayResponse_cell(sqrtNt, Ncl, Nray, std);
        Arcell = ArrayResponse_cell(sqrtNr, Ncl, Nray, std);
        
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
        I1(l,i) = SS(1, NtRF, NrRF, H, At, Ar, sigma_n, rho(i));
        I2(l,i) = SS(2, NtRF, NrRF, H, At, Ar, sigma_n, rho(i));
        
        % unconstrained - Unconstrained Precoding / Decoding
        I3(l,i) = unconstrained (1, H, sigma_n, rho(i));
        I4(l,i) = unconstrained (2, H, sigma_n, rho(i));    
    end
end

I1mean = mean(I1);
I2mean = mean(I2);
I3mean = mean(I3);
I4mean = mean(I4);

x = plot(rhodB, I1mean, rhodB, I2mean, rhodB, I3mean, rhodB, I4mean);

x(1).LineWidth = 2;
x(2).LineWidth = 2;
x(3).LineWidth = 2;
x(4).LineWidth = 2;

x(1).Marker = 'o';
x(2).Marker = 'o';
x(3).Marker = 's';
x(4).Marker = 's';

legend('SS Precoding & Combining, N_{s}=1', 'SS Precoding & Combining, N_{s}=2', 'Optimal Unstrained Precoding, N_{s}=1', 'Optimal Unstrained Precoding, N_{s}=2', 'Location', 'northwest');
title('Spatially Sparse Precoding  vs Unconstrained Optimum Precoding');
xlabel('SNR (dB)');
ylabel('Spectral Efficiency (bits/s/Hz)');
grid on;

end