function [SNR1, SNR2] = SNR_QoS (sigma, H, D)

% QoS Based Design %

% System Size Def
B = 2;
[Mr, ~] = size(H);

% Noise Matrix Def
Rnn = sigma.^2 * eye(Mr);

[~, LAMBDA] = eig(H' * Rnn^(-1) * H);

LAMBDA = sort(diag(LAMBDA), 'descend');

LAMBDA0 = zeros(B,B);
for i = 1 : B
    LAMBDA0(i,i) = LAMBDA(i);
end

% QoS Optimization
gamma = 1 / trace(D * LAMBDA0^(-1));            % p0 = 1
PHIf = gamma^0.5 * D^0.5 * LAMBDA0^(-0.5);

% SNR, Bit Rate
SNR = PHIf^2 * LAMBDA0;
SNR1 = SNR(1,1);
SNR2 = SNR(2,2);

end