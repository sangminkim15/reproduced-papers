function [objective_value] = algorithm_MCEU(p, H)

%% 1. Initialization \Wv_n and \mathbf{\Phi}
RAND11 = normrnd(0,1/sqrt(2),[p.N_t,p.d]) + 1i * normrnd(0,1/sqrt(2),[p.N_t,p.d]);
RAND12 = normrnd(0,1/sqrt(2),[p.N_t,p.d]) + 1i * normrnd(0,1/sqrt(2),[p.N_t,p.d]);
RAND13 = normrnd(0,1/sqrt(2),[p.N_t,p.d]) + 1i * normrnd(0,1/sqrt(2),[p.N_t,p.d]);
RAND21 = normrnd(0,1/sqrt(2),[p.N_t,p.d]) + 1i * normrnd(0,1/sqrt(2),[p.N_t,p.d]);
RAND22 = normrnd(0,1/sqrt(2),[p.N_t,p.d]) + 1i * normrnd(0,1/sqrt(2),[p.N_t,p.d]);
RAND23 = normrnd(0,1/sqrt(2),[p.N_t,p.d]) + 1i * normrnd(0,1/sqrt(2),[p.N_t,p.d]);
RAND31 = normrnd(0,1/sqrt(2),[p.N_t,p.d]) + 1i * normrnd(0,1/sqrt(2),[p.N_t,p.d]);
RAND32 = normrnd(0,1/sqrt(2),[p.N_t,p.d]) + 1i * normrnd(0,1/sqrt(2),[p.N_t,p.d]);
RAND33 = normrnd(0,1/sqrt(2),[p.N_t,p.d]) + 1i * normrnd(0,1/sqrt(2),[p.N_t,p.d]);

Wprev.W11 = sqrt(p.P_max / 3) * RAND11;         Wprev.W12 = sqrt(p.P_max / 3) * RAND12;         Wprev.W13 = sqrt(p.P_max / 3) * RAND13;
Wprev.W21 = sqrt(p.P_max / 3) * RAND21;         Wprev.W22 = sqrt(p.P_max / 3) * RAND22;         Wprev.W23 = sqrt(p.P_max / 3) * RAND23;
Wprev.W31 = sqrt(p.P_max / 3) * RAND31;         Wprev.W32 = sqrt(p.P_max / 3) * RAND32;         Wprev.W33 = sqrt(p.P_max / 3) * RAND33;

PHIprev = eye(p.M, p.M);

U = Uopt_MCEU (p, H, Wprev, PHIprev);

Q = Qopt_MCEU (p, H, Wprev, PHIprev);

[W, Rtemp] = Wopt_MCEU (p, H, U, Q, PHIprev);

[PHI, R] = PHIopt_MCEU (p, H, U, Q, W);

if real(Rtemp) > real(R)
    PHI = PHIprev;
end

while true
    Wprev = W; PHIprev = PHI;
    Rprev = R;
    
    %% 2. Calculte \Wv^{opt} from (15)
    U = Uopt_MCEU (p, H, Wprev, PHIprev);
    
    %% 3. Calculte \Qv^{opt} from (16)
    Q = Qopt_MCEU (p, H, Wprev, PHIprev);
    
    %% 4. Calculate \Wv^{opt}_n from Algorithm 1.
    [W, Rtemp] = Wopt_MCEU (p, H, U, Q, PHIprev);
    
    %% 5. Calculate \Phi^{opt}_n from Algorithm 2.
    [PHI, R] = PHIopt_MCEU (p, H, U, Q, W);
    
    if real(Rtemp) > real(R)
        PHI = PHIprev;
    end
    
    %% Stopping criteria
    if abs(R - Rprev) < p.eta
        break;
    end
end

objective_value = R;

end

