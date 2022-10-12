function [WSR] = TxZF (u, H, Etx)

% I/O
% H     Channel Matrix, size QK * P
% u     Weight, size 1 * K
% Etx   Power Constraint

% WSR   Weight Sum Rate

[K, ~] = size(H);

Hmatrix = cell2mat(H);

Bmatrix = Hmatrix' * (Hmatrix*Hmatrix')^(-1);
temp = norm(Bmatrix, 'fro').^2;

Bmatrix = (Etx / temp).^0.5 * Bmatrix;

B = matrixdiv(Bmatrix, K);

E = cell(1, K);
R = zeros(1, K);

for i = 1 : K
    [E{1,i}, ~] = WMMSE(u(i), H{i,1}, B, i);
    R(i) = log(det(E{1,i}^(-1))) / log(2);
end

WSR = sum(u.*R);

end