function [PHI] = algorithm_2 (p, H, U, Q, W1, W2, PHIprev)

HW = H.bs1_ue1 * W1 + H.bs2_ue1 * W2;
GW = H.bs1_IRS * W1 + H.bs2_IRS * W2;

A = H.IRS_ue1' * U * Q * U' * H.IRS_ue1;        Etilde = GW * GW';
B = H.IRS_ue1' * U * Q * U' * GW';              D = H.IRS_ue1' * U * Q * U' * HW * GW';

AE = A .* transpose(Etilde);
eigmax = max(eig(AE));

phiprev = diag(PHIprev);

z = diag(D - B);
q = z - (eigmax * eye(size(AE)) - AE) * phiprev;
phi = - exp(1i .* angle(q));

S = eigmax * (phi' * phi) - 2 * real(phi' * (eigmax * eye(size(AE)) - AE) * phiprev) + phiprev' * (eigmax * eye(size(AE)) - AE) * phiprev + z' * phi + phi' * z;

while true
    phiprev = phi;
    Sprev = S;
    
    q = z - (eigmax * eye(size(AE)) - AE) * phiprev;
    phi = -exp(1i .* angle(q));
    
    S = eigmax * (phi' * phi) - 2 * real(phi' * (eigmax * eye(size(AE)) - AE) * phiprev) + phiprev' * (eigmax * eye(size(AE)) - AE) * phiprev + z' * phi + phi' * z;
    
    if abs(S - Sprev) < p.zeta
        break
    end
end

PHI = zeros(length(phi), length(phi));
for i = 1 : length(phi)
    PHI(i,i) = phi(i);
end

end