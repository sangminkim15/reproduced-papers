function [Ek, Wk] = WMMSE (uk, Hk, B, k)

% MMSE Weight Matrix %

% I/O
% Hk    Channel Matrix of User k, size Q * P
% B     Beamforming Matrix, size P * QK
%   B = [B1 B2 ... Bk ... BK]
% k     User No.    

% Wk    MMSE Weight Matrix of User k

[Q, P] = size(Hk);
[~, K] = size(B);

Bk = zeros(P,Q);
Y = B{1,k};
for i = 1 : P
    for j = 1 : Q
        Bk(i,j) = Y(i,j);
    end
end

Rvkvk = eye(Q,Q) - Hk * (Bk * Bk') * Hk';

for l = 1 : K
    Rvkvk = Rvkvk + Hk * (B{1,l} * B{1,l}') * Hk';
end

Ek = (eye(Q,Q) + Bk' * Hk' * Rvkvk^(-1) * Hk * Bk)^(-1);
Wk = uk * Ek^(-1);

end