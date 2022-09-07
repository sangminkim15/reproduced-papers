function [CELL] = matrixdiv (MATRIX, K)

% Divide P * QK matrix to K P * Q size matrices

CELL = cell(1,K);

[P, QK] = size(MATRIX);
Q = QK / K;

for i = 1 : K
    X = zeros(P, Q);
    for j = 1 : P
        for l = 1 : Q
            X(j,l) = MATRIX(j, l + Q * (i-1));
        end
    end
    CELL{1,i} = X;
end

end