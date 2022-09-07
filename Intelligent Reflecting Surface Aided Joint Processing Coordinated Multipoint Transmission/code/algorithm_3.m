function [objective_value] = algorithm_3(p,H)

warning off

%% 1. Initialization \Wv_n and \mathbf{\Phi}
RAND1 = normrnd(0,1/sqrt(2),[p.N_t,p.d]) + 1i * normrnd(0,1/sqrt(2),[p.N_t,p.d]);
RAND2 = normrnd(0,1/sqrt(2),[p.N_t,p.d]) + 1i * normrnd(0,1/sqrt(2),[p.N_t,p.d]);
W1prev = (sqrt(p.P_max)/norm(RAND1, 'fro')) * RAND1;       W2prev = (sqrt(p.P_max)/norm(RAND2, 'fro')) * RAND2;

PHIprev = eye(p.M, p.M);

U = Uopt (p, H, W1prev, W2prev, PHIprev);

Q = Qopt (p, H, W1prev, W2prev, PHIprev);

W = algorithm_1 (p, H, U, Q, PHIprev);

sz = size(W);
W1 = W(1:sz(1)/2, 1:sz(2));
W2 = W(sz(1)/2+1:sz(1), 1:sz(2));

PHI = algorithm_2 (p, H, U, Q, W1, W2, PHIprev);

R = real(log(det(Q)));

while true
    Wprev = W; PHIprev = PHI;
    
    szprev = size(Wprev);
    
    W1prev = Wprev(1:szprev(1)/2, 1:szprev(2));
    W2prev = Wprev(szprev(1)/2+1:szprev(1), 1:szprev(2));
    
    %% 2. Calculte \Wv^{opt} from (15)
    U = Uopt (p, H, W1prev, W2prev, PHIprev);
    
    %% 3. Calculte \Qv^{opt} from (16)
    Q = Qopt (p, H, W1prev, W2prev, PHIprev);
    
    %% 4. Calculate \Wv^{opt}_n from Algorithm 1.
    W = algorithm_1 (p, H, U, Q, PHIprev);
    
    sz = size(W);
    W1 = W(1:sz(1)/2, 1:sz(2));
    W2 = W(sz(1)/2+1:sz(1), 1:sz(2));
    
    %% 5. Calculate \Phi^{opt}_n from Algorithm 2.
    PHI = algorithm_2 (p, H, U, Q, W1, W2, PHIprev);
    
    %% Stopping criteria
    Rprev = R;          R = real(log(det(Q)));
    
    if abs(R - Rprev) < p.eta
        break;
    end
end

objective_value = R;

end