function [next] = stsf(curr, lambda)

% Soft Thresholding Shrinkage Function (STSF)

[M, ~] = size(curr);
next = zeros(M, 1);

for i = 1 : M
    next(i) = sign(curr(i)) * max(abs(curr(i)) - lambda, 0);
end

end