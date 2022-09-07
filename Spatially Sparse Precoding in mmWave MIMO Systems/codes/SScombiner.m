function [WRF, WBB] = SScombiner (Ns, NrRF, H, Ar, FRF, FBB, sigma_n, rho)

% Spatially Sparse MMSE Combining (SS Combining)

% I/O
% Nr    # of receiver antennas
% NrRF  # of receiver RF chains
% WRF   RF Combiner (Nr * NrRF)
% WBB   Baseband Combiner (NrRF * Ns)

WRFnext = zeros(0, 0);
Wmmse = WMMSE(Ns, H, FRF, FBB, sigma_n, rho);

Nr = size(H * FRF * (FBB * FBB') * FRF' * H');
Eyy = (rho/Ns) * H * FRF * (FBB * FBB') * FRF' * H' + sigma_n.^2 * eye(Nr);

Wres = Wmmse;

for i = 1 : NrRF
    WRFprev = WRFnext;
    PSI = Ar' * Eyy * Wres;
    k = diag(PSI * PSI') == max(diag(PSI * PSI'));
    
    WRFnext = [WRFprev Ar(:,k)];
    WBB = (WRFnext' * Eyy * WRFnext) \ (WRFnext' * Eyy * Wmmse);
    
    Wres = (Wmmse - WRFnext * WBB) / norm(Wmmse - WRFnext * WBB, 'fro');
end

WRF = WRFnext;

end