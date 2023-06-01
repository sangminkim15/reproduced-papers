function [X] = fista(A, y, T, lambda, beta)

% Fast ISTA (FISTA)

[~, N] = size(A);
i = 0;
curr = zeros(N,1);
prev = curr;
X = zeros(N,T);

while i < T
    v = y - A * curr;
    next = stsf(curr + beta * transpose(A) * v + (i-2)/(i+1) * (curr - prev), lambda);
    prev = curr;
    curr = next;
    i = i + 1;
    X(:,i) = curr;
end

end