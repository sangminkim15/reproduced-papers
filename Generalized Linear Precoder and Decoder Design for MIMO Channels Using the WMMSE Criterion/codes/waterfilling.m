function [LAMBDA1, PHIf, SNR, R, M] = waterfilling ()

% Waterfilling Method : W = LAMBDA %

% System Size Def
B = 5;
Mt = 5;
Mr = 5;

% Channel / Noise Matrix Def
H = sqrt(1/2) * (randn(Mr,Mt) + 1i * randn(Mr,Mt));
sigma = 1e-1;
Rnn = sigma.^2 * eye(Mr);

% Eigenvalue Decomposition / SOrting
[~, LAMBDA] = eig(H' * Rnn^(-1) * H);

LAMBDA = sort(diag(LAMBDA), 'descend');

LAMBDA0 = zeros(B,B);
for i = 1 : B
    LAMBDA0(i,i) = LAMBDA(i);
end

% Precoder / Decoder Optimization
[~, LAMBDA1, PHIf, ~] = optimum(B, LAMBDA0, LAMBDA0, 1);

% SNR, Bit Rate
SNR = PHIf^2 * LAMBDA1;
R = log(ones(B,B) + SNR) / log(2);
M = 2.^R - (ones(B,B) - eye(B));

end
