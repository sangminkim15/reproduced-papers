function [X] = ista(A, y, T, lambda, beta)

% Iterative Shrinkage/Thresholding Algorithm (ISTA)

[~, N] = size(A);
i = 0;
x = zeros(N,1);
X = zeros(N,T);

while i < T
    v = y - A * x;
    x = stsf(x + beta * transpose(A) * v, lambda);
    i = i + 1;
    X(:,i) = x;
end

end