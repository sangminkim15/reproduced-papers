function [I] = spectraleff (Ns, H, FRF, FBB, WRF, WBB, sigma_n, rho)

% Spectral Effieciency

Rn = sigma_n.^2 * WBB' * (WRF' * WRF) * WBB;
I = log(det(eye(Ns) + (rho/Ns) * (Rn \ WBB' * WRF' * H * FRF * (FBB * FBB') * FRF' * H' * WRF * WBB)))/log(2);

end