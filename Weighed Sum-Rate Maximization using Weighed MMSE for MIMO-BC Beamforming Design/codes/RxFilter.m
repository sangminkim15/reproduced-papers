function [Ak] = RxFilter (Hk, B, k)

% MMSE Receiver Filter %

% I/O
% Hk    Channel Matrix of User k, size Q * P
% B     Beamforming Matrix, size P * QK
%   B = [B1 B2 ... Bk ... BK]
% k     User No.    

% Ak    MMSE Receiver Filter of User k

[Q, P] = size(Hk);
[~, K] = size(B);

Bk = zeros(P,Q);
Y = B{1,k};
for i = 1 : P
    for j = 1 : Q
        Bk(i,j) = Y(i,j);
    end
end

X = eye(Q,Q);

for l = 1 : K
    X = X + Hk * (B{1,l} * B{1,l}') * Hk';      % X = H * (Bk * Bk') * H' + Rvkvk
end

Ak = Bk' * Hk' * X^(-1);

end