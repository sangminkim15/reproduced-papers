function [WSR] = TxMF (u, H, Etx)

% I/O
% H     Channel Matrix, size QK * P
% u     Weight, size 1 * K
% Etx   Power Constraint

% WSR   Weight Sum Rate

% Initialization
[K, ~] = size(H);

B = cell(1,K);
tr = 0;
for i = 1 : K
    B{1,i} = H{i,1}';
    tr = tr + trace(B{1,i} * B{1,i}');
end

for i = 1 : K
    B{1,i} = (Etx / tr).^0.5 * B{1,i};
end

E = cell(1, K);
R = zeros(1, K);

for i = 1 : K
    [E{1,i}, ~] = WMMSE(u(i), H{i,1}, B, i);
    R(i) = log(det(E{1,i}^(-1))) / log(2);
end

WSR = sum(u.*R);

end