function [W, Rate] = Wopt_MCEU (p, H, U, Q, PHI)

H11 = H.bs1_ue1 + H.IRS_ue1 * PHI * H.bs1_IRS;      H21 = H.bs2_ue1 + H.IRS_ue1 * PHI * H.bs2_IRS;      H31 = H.bs3_ue1 + H.IRS_ue1 * PHI * H.bs3_IRS;
H12 = H.bs1_ue2 + H.IRS_ue1 * PHI * H.bs2_IRS;      H22 = H.bs2_ue2 + H.IRS_ue2 * PHI * H.bs2_IRS;      H32 = H.bs3_ue2 + H.IRS_ue2 * PHI * H.bs3_IRS;
H13 = H.bs1_ue3 + H.IRS_ue1 * PHI * H.bs3_IRS;      H23 = H.bs2_ue3 + H.IRS_ue3 * PHI * H.bs2_IRS;      H33 = H.bs3_ue3 + H.IRS_ue3 * PHI * H.bs3_IRS;

H1 = [H11 H21 H31];                                 H2 = [H12 H22 H32];                                 H3 = [H13 H23 H33];

cvx_begin

variable W11(p.N_t, p.d) complex
variable W12(p.N_t, p.d) complex
variable W13(p.N_t, p.d) complex
variable W21(p.N_t, p.d) complex
variable W22(p.N_t, p.d) complex
variable W23(p.N_t, p.d) complex
variable W31(p.N_t, p.d) complex
variable W32(p.N_t, p.d) complex
variable W33(p.N_t, p.d) complex
variable R

W1 = [W11 ; W21 ; W31];         W2 = [W12 ; W22 ; W32];         W3 = [W13 ; W23 ; W33];

eta1 = [W11(:) ; W12(:) ; W13(:)];
eta2 = [W21(:) ; W22(:) ; W23(:)];
eta3 = [W31(:) ; W32(:) ; W33(:)];

omega1 = [vec(W1'*H1'*U.U1*Q.Q1^0.5 - Q.Q1^0.5) ; vec(W2'*H1'*U.U1*Q.Q1^0.5) ; vec(W3'*H1'*U.U1*Q.Q1^0.5)];
omega2 = [vec(W1'*H2'*U.U2*Q.Q2^0.5) ; vec(W2'*H2'*U.U2*Q.Q2^0.5 - Q.Q2^0.5) ; vec(W3'*H2'*U.U2*Q.Q2^0.5)];
omega3 = [vec(W1'*H3'*U.U3*Q.Q3^0.5) ; vec(W2'*H3'*U.U3*Q.Q3^0.5) ; vec(W3'*H3'*U.U3*Q.Q3^0.5 - Q.Q3^0.5)];

maximize (R)

subject to
norm(eta1) <= p.P_max;
norm(eta2) <= p.P_max;
norm(eta3) <= p.P_max;
norm(omega1) <= sqrt(real(log(det(Q.Q1)) + p.d - R - p.np * trace(Q.Q1 * U.U1' * U.U1)));
norm(omega2) <= sqrt(real(log(det(Q.Q2)) + p.d - R - p.np * trace(Q.Q2 * U.U2' * U.U2)));
norm(omega3) <= sqrt(real(log(det(Q.Q3)) + p.d - R - p.np * trace(Q.Q3 * U.U3' * U.U3)));

cvx_end

W.W11 = W11;        W.W12 = W12;        W.W13 = W13;
W.W21 = W21;        W.W22 = W22;        W.W23 = W23;
W.W31 = W31;        W.W32 = W32;        W.W33 = W33;

Rate = R;

end