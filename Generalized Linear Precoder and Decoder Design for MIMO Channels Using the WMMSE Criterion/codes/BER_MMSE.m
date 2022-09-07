function [BER, avgBER] = BER_MMSE (sigma, H, B)

% MMSE Bit Error Rate : W = I %

% System Size Def
[Mr, ~] = size(H);

% Noise Def
Rnn = sigma.^2 * eye(Mr);

% Eigenvalue Decomposition / Sorting
[~, LAMBDA] = eig(H' * Rnn^(-1) * H);

LAMBDA = sort(diag(LAMBDA), 'descend');

LAMBDA0 = zeros(B,B);
for i = 1 : B
    LAMBDA0(i,i) = LAMBDA(i);
end

% Precoder / Decoder Optimization
[~, LAMBDA1, PHIf, ~] = optimum(B, LAMBDA0, eye(B), 1);

% SNR, BER
SNR = PHIf^2 * LAMBDA1;
BER = 2 * (0.5 - (0.5) * erf(sqrt(real(SNR))));
avgBER = mean(diag(BER), 'all');

end