function [B] = Beamforming (H, A, W, Etx)

% Beamforming Matrix %

% I/O
% H     Channel Matrix, size QK * P
%   H = [H1 ; H2 ; ... ; Hk ; ... ; HK]
% A     MMSE Receive Filter Matrix, size QK * QK
% W     MMSE Weight Matrix, size QK * QK
% Etx   Power Constraint    

% B     Beamforming Matrix

[~, P] = size(H);

B = (H' * A' * W * A * H + (trace(W * (A * A')) / Etx) * eye(P, P))^(-1) * H' * A' * W;
b = (Etx / trace(B * B'))^(0.5);
B = b * B;

end