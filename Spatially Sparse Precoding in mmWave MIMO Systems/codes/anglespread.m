function [Issmean, Iunmean] = anglespread (Ns, Nt, Nr, NtRF, NrRF, Ncl, Nray, std)

% # of antennas
sqrtNt = sqrt(Nt);
sqrtNr = sqrt(Nr);

% SNR fixed
sigma_n = 1;    % 0dB
rho = 1;

Iss = zeros(1000, length(std));     % SS (Spatially Sparse)
Iun = zeros(1000, length(std));     % Unconstrained

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
        Iss(l,i) = SS(Ns, NtRF, NrRF, H, At, Ar, sigma_n, rho);
        
        % unconstrained - Unconstrained Precoding / Decoding
        Iun(l,i) = unconstrained(Ns, H, sigma_n, rho);
    end
end

Issmean = mean(Iss);
Iunmean = mean(Iun);

end