function [Wmmse] = WMMSE(Ns, H, FRF, FBB, sigma_n, rho)

% Calculate Wmmse for SS combining

Esy = (sqrt(rho)/Ns) * FBB' * FRF' * H';

[Nr, ~] = size(H * FRF * (FBB * FBB') * FRF' * H');
Eyy = (rho/Ns) * H * FRF * (FBB * FBB') * FRF' * H' + sigma_n.^2 * eye(Nr);

Wmmse = (Esy / Eyy)';

end