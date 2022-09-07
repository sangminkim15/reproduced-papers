function [Q] = Qopt_MCEU (p, H, W, PHI)

H11 = H.bs1_ue1 + H.IRS_ue1 * PHI * H.bs1_IRS;      H21 = H.bs2_ue1 + H.IRS_ue1 * PHI * H.bs2_IRS;      H31 = H.bs3_ue1 + H.IRS_ue1 * PHI * H.bs3_IRS;
H12 = H.bs1_ue2 + H.IRS_ue1 * PHI * H.bs2_IRS;      H22 = H.bs2_ue2 + H.IRS_ue2 * PHI * H.bs2_IRS;      H32 = H.bs3_ue2 + H.IRS_ue2 * PHI * H.bs3_IRS;
H13 = H.bs1_ue3 + H.IRS_ue1 * PHI * H.bs3_IRS;      H23 = H.bs2_ue3 + H.IRS_ue3 * PHI * H.bs2_IRS;      H33 = H.bs3_ue3 + H.IRS_ue3 * PHI * H.bs3_IRS;

H1 = [H11 H21 H31];                                 H2 = [H12 H22 H32];                                 H3 = [H13 H23 H33];

W1 = [W.W11 ; W.W21 ; W.W31];                       W2 = [W.W12 ; W.W22 ; W.W32];                       W3 = [W.W13 ; W.W23 ; W.W33];

J1 = (H1 * (W1 * W1' + W2 * W2' + W3 * W3') * H1') + p.np * eye(p.N_t);
J2 = (H2 * (W1 * W1' + W2 * W2' + W3 * W3') * H2') + p.np * eye(p.N_t);
J3 = (H3 * (W1 * W1' + W2 * W2' + W3 * W3') * H3') + p.np * eye(p.N_t);

Q.Q1 = inv(eye(p.d) - ((W1' * H1') / J1) * H1 * W1);
Q.Q2 = inv(eye(p.d) - ((W2' * H2') / J2) * H2 * W2);
Q.Q3 = inv(eye(p.d) - ((W3' * H3') / J3) * H3 * W3);

end