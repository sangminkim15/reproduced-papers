function [Issmean, Iunmean] = SNR (Ns, Nt, Nr, NtRF, NrRF, Ncl, Nray, rhodB)

% # of antennas
sqrtNt = sqrt(Nt);
sqrtNr = sqrt(Nr);

% Angle Spread Fixed
std = pi/24;

% SNR
sigma_n = 1;                       % 0dB
rho = 10.^(rhodB/10);

Iss = zeros(1000, length(rho));     % SS (Spatially Sparse)
Iun = zeros(1000, length(rho));     % Unconstrained

for i = 1 : length(rho)
    for l = 1 : 1000
        Atcell =ArrayResponse_cell(sqrtNt, Ncl, Nray, std);
        Arcell =ArrayResponse_cell(sqrtNr, Ncl, Nray, std);
        
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
        Iss(l,i) = SS(Ns, NtRF, NrRF, H, At, Ar, sigma_n, rho(i));
        
        % unconstrained - Unconstrained Precoding / Decoding
        Iun(l,i) = unconstrained(Ns, H, sigma_n, rho(i));
    end
end

Issmean = mean(Iss);
Iunmean = mean(Iun);

end