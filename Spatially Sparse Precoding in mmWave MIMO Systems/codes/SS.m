function [I] = SS(Ns, NtRF, NrRF, H, At, Ar, sigma_n, rho)

% Spatially Sparse 

[FRF, FBB] = SSprecoder(Ns, NtRF, H, At);
[WRF, WBB] = SScombiner(Ns, NrRF, H, Ar, FRF, FBB, sigma_n, rho);

I = spectraleff(Ns, H, FRF, FBB, WRF, WBB, sigma_n, rho);

end