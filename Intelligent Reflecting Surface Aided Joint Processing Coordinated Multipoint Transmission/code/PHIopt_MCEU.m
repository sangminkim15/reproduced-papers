function [PHI, Rate] = PHIopt_MCEU (p, H, U, Q, W)

L1.L1 = H.bs1_IRS * W.W11 + H.bs2_IRS * W.W21 + H.bs3_IRS * W.W31;
L1.L2 = H.bs1_IRS * W.W12 + H.bs2_IRS * W.W22 + H.bs3_IRS * W.W32;
L1.L3 = H.bs1_IRS * W.W13 + H.bs2_IRS * W.W23 + H.bs3_IRS * W.W33;

L2.L11 = H.bs1_ue1 * W.W11 + H.bs2_ue1 * W.W21 + H.bs3_ue1 * W.W31;
L2.L12 = H.bs1_ue1 * W.W12 + H.bs2_ue1 * W.W22 + H.bs3_ue1 * W.W32;
L2.L13 = H.bs1_ue1 * W.W13 + H.bs2_ue1 * W.W23 + H.bs3_ue1 * W.W33;

L2.L21 = H.bs1_ue2 * W.W11 + H.bs2_ue2 * W.W21 + H.bs3_ue2 * W.W31;
L2.L22 = H.bs1_ue2 * W.W12 + H.bs2_ue2 * W.W22 + H.bs3_ue2 * W.W32;
L2.L23 = H.bs1_ue2 * W.W13 + H.bs2_ue2 * W.W23 + H.bs3_ue2 * W.W33;

L2.L31 = H.bs1_ue3 * W.W11 + H.bs2_ue3 * W.W21 + H.bs3_ue3 * W.W31;
L2.L32 = H.bs1_ue3 * W.W12 + H.bs2_ue3 * W.W22 + H.bs3_ue3 * W.W32;
L2.L33 = H.bs1_ue3 * W.W13 + H.bs2_ue3 * W.W23 + H.bs3_ue3 * W.W33;

A.A1 = H.IRS_ue1' * U.U1 * Q.Q1 * U.U1' * H.IRS_ue1;
A.A2 = H.IRS_ue2' * U.U2 * Q.Q2 * U.U2' * H.IRS_ue2;
A.A3 = H.IRS_ue3' * U.U3 * Q.Q3 * U.U3' * H.IRS_ue3;

B.B1 = H.IRS_ue1' * U.U1 * Q.Q1 * L1.L1';
B.B2 = H.IRS_ue2' * U.U2 * Q.Q2 * L1.L2';
B.B3 = H.IRS_ue3' * U.U3 * Q.Q3 * L1.L3';

E = L1.L1 * L1.L1' + L1.L2 * L1.L2' + L1.L3 * L1.L3';

D.D1 = H.IRS_ue1' * U.U1 * Q.Q1 * U.U1' * (L2.L11 * L1.L1' + L2.L12 * L1.L2' + L2.L13 * L1.L3');
D.D2 = H.IRS_ue2' * U.U2 * Q.Q2 * U.U2' * (L2.L21 * L1.L1' + L2.L22 * L1.L2' + L2.L23 * L1.L3');
D.D3 = H.IRS_ue3' * U.U3 * Q.Q3 * U.U3' * (L2.L31 * L1.L1' + L2.L32 * L1.L2' + L2.L33 * L1.L3');

c1.c1 = trace(Q.Q1 * L2.L11' * U.U1);
c1.c2 = trace(Q.Q2 * L2.L22' * U.U2);
c1.c3 = trace(Q.Q3 * L2.L33' * U.U3);

c2.c1 = trace((L2.L11 * L2.L11' + L2.L12 * L2.L12' + L2.L13 * L2.L13') * U.U1 * Q.Q1 * U.U1');
c2.c2 = trace((L2.L21 * L2.L21' + L2.L22 * L2.L22' + L2.L23 * L2.L23') * U.U2 * Q.Q2 * U.U2');
c2.c3 = trace((L2.L31 * L2.L31' + L2.L32 * L2.L32' + L2.L33 * L2.L33') * U.U1 * Q.Q1 * U.U1');

const.c1 = real(log(det(Q.Q1))) + p.d + 2 * real(c1.c1) - c2.c1 - trace(Q.Q1 * (p.variance .* U.U1' * U.U1 + eye(p.d, p.d)));
const.c2 = real(log(det(Q.Q2))) + p.d + 2 * real(c1.c2) - c2.c2 - trace(Q.Q2 * (p.variance .* U.U2' * U.U2 + eye(p.d, p.d)));
const.c3 = real(log(det(Q.Q3))) + p.d + 2 * real(c1.c3) - c2.c3 - trace(Q.Q3 * (p.variance .* U.U3' * U.U3 + eye(p.d, p.d)));

z.z1 = diag (D.D1 - B.B1);
z.z2 = diag (D.D2 - B.B2);
z.z3 = diag (D.D3 - B.B3);

PSI.PSI1 = [A.A1.*transpose(E) z.z1 ; z.z1' 0];
PSI.PSI2 = [A.A2.*transpose(E) z.z2 ; z.z2' 0];
PSI.PSI3 = [A.A3.*transpose(E) z.z3 ; z.z3' 0];

THETA = zeros(p.M+1, p.M+1);

cvx_begin sdp
variable THETA(p.M + 1, p.M + 1) semidefinite
variable R

maximize (R)

subject to
real(trace(PSI.PSI1 * THETA)) <= real(const.c1 - R);
real(trace(PSI.PSI2 * THETA)) <= real(const.c2 - R);
real(trace(PSI.PSI3 * THETA)) <= real(const.c3 - R);
diag(THETA) == ones(p.M+1,1);

cvx_end

[U, SIGMA] = eigs(THETA);
sz = size(SIGMA);

phibararray = zeros(0, 0);
Rarray = zeros(0, 0);

for idx = 1 : 1000
    randvec = (1 / sqrt(2)) * (randn(sz(2), 1) + 1i * randn(sz(2), 1));
    phitilde = U * SIGMA^0.5 * randvec;

    phibar = exp(1i .* angle(phitilde ./ phitilde(p.M + 1)));
    phibararray = horzcat(phibararray, phibar);
    
    PHIbar = zeros(p.M, p.M);
    for jdx = 1 : p.M
        PHIbar(jdx, jdx) = phibar(jdx);
    end
    
    H11 = H.bs1_ue1 + H.IRS_ue1 * PHIbar * H.bs1_IRS;      H21 = H.bs2_ue1 + H.IRS_ue1 * PHIbar * H.bs2_IRS;      H31 = H.bs3_ue1 + H.IRS_ue1 * PHIbar * H.bs3_IRS;
    H12 = H.bs1_ue2 + H.IRS_ue1 * PHIbar * H.bs2_IRS;      H22 = H.bs2_ue2 + H.IRS_ue2 * PHIbar * H.bs2_IRS;      H32 = H.bs3_ue2 + H.IRS_ue2 * PHIbar * H.bs3_IRS;
    H13 = H.bs1_ue3 + H.IRS_ue1 * PHIbar * H.bs3_IRS;      H23 = H.bs2_ue3 + H.IRS_ue3 * PHIbar * H.bs2_IRS;      H33 = H.bs3_ue3 + H.IRS_ue3 * PHIbar * H.bs3_IRS;

    H1 = [H11 H21 H31];                                 H2 = [H12 H22 H32];                                 H3 = [H13 H23 H33];

    W1 = [W.W11 ; W.W21 ; W.W31];                       W2 = [W.W12 ; W.W22 ; W.W32];                       W3 = [W.W13 ; W.W23 ; W.W33];
    
    J1 = (H1 * (W1 * W1' + W2 * W2' + W3 * W3') * H1') + p.np * eye(p.N_t);
    J2 = (H2 * (W1 * W1' + W2 * W2' + W3 * W3') * H2') + p.np * eye(p.N_t);
    J3 = (H3 * (W1 * W1' + W2 * W2' + W3 * W3') * H3') + p.np * eye(p.N_t);

    Q1 = inv(eye(p.d) - ((W1' * H1') / J1) * H1 * W1);
    Q2 = inv(eye(p.d) - ((W2' * H2') / J2) * H2 * W2);
    Q3 = inv(eye(p.d) - ((W3' * H3') / J3) * H3 * W3);
    
    R = min([log(det(Q1)), log(det(Q2)), log(det(Q3))]);
    
    Rarray = horzcat(Rarray, R);
end

maxindex = find(Rarray == max(Rarray));

phiopt = phibararray(:, maxindex);

Rate = Rarray(maxindex);

PHI = zeros(p.M, p.M);

for kdx = 1 : p.M
    PHI(kdx, kdx) = phiopt(kdx);
end

end