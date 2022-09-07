function [WSR] = WSRBF_WMMSE2 (u, H, Etx)

% Algorithm 2
%   Simple Transmit Matched Filter Design & 10 Iterations

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

[~,Q] = size(B{1,1});

A = cell(K,K);
W = cell(K,K);

% 10 Iterations
n = 0;
while n < 10
    for j = 1 : K
        for l = 1 : K
            if j == l
                A{j,l} = RxFilter(H{j,1}, B, j);
                [~, W{j,l}] = WMMSE(u(j), H{j,1}, B, j);
            else
                A{j,l} = zeros(Q,Q);
                W{j,l} = zeros(Q,Q);
            end
        end
    end
    
    Hmatrix = cell2mat(H);
    Amatrix = cell2mat(A);
    Wmatrix = cell2mat(W);
    
    Bmatrix = Beamforming(Hmatrix, Amatrix, Wmatrix, Etx);
    
    B = matrixdiv(Bmatrix, K);
    
    n = n + 1;
end

E = cell(1, K);
R = zeros(1, K);

for i = 1 : K
    [E{1,i}, ~] = WMMSE(u(i), H{i,1}, B, i);
    R(i) = log(det(E{1,i}^(-1))) / log(2);
end

WSR = sum(u.*R);

end