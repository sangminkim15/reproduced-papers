function [I] = unconstrained (Ns, H, sigma_n, rho)

FrfFbb = Fopt(Ns, H);
Wmmse = Wopt(Ns, H, FrfFbb, sigma_n, rho);

Rn = sigma_n.^2 * (Wmmse' * Wmmse);
I = log(det(eye(Ns) + (rho/Ns) * (Rn \ Wmmse' * H * (FrfFbb * FrfFbb') * H' * Wmmse)))/log(2);

end