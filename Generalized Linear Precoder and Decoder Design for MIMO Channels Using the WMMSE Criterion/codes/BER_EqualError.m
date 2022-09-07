function [avgBER] = BER_EqualError (sigma, H, B)

% Equal-Error Design %

% System Size Def
[Mr, ~] = size(H);

% Noise Matrix Def
Rnn = sigma.^2 * eye(Mr);

[~, LAMBDA] = eig(H' * Rnn^(-1) * H);

LAMBDA = sort(diag(LAMBDA), 'descend');

LAMBDA0 = zeros(B,B);
for i = 1 : B
    LAMBDA0(i,i) = LAMBDA(i);
end

% Equal-Error Optimization
D = (1/B) * eye(B);
gamma = 1 / trace(D * LAMBDA0^(-1));            % p0 = 1
PHIf = gamma^0.5 * D^0.5 * LAMBDA0^(-0.5);

% SNR, Bit Rate
SNR = PHIf^2 * LAMBDA0;
BER = 2 * (0.5 - (0.5) * erf(sqrt(real(SNR))));
avgBER = mean(diag(BER), 'all');

end