function [L] = Lagrangian (p, H, U, Q, W1, W2, PHI, mu)

H1b = H.bs1_ue1 + H.IRS_ue1 * PHI * H.bs1_IRS;      H2b = H.bs2_ue1 + H.IRS_ue1 * PHI * H.bs2_IRS;
HW = H1b * W1 + H2b * W2;

L = trace(HW * HW' * U * Q * U') - trace(Q * U' * HW) - trace(Q * HW' * U) + mu(1) * (norm(W1, 'fro') - p.P_max) + mu(2) * (norm(W2, 'fro') - p.P_max);

end