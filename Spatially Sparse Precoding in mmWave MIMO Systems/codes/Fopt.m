function [V1] = Fopt (Ns, H)

% Calculates Unconstrained Unitary Precoder

[~, Nt] = size(H);
V1 = zeros(Nt, Ns);

[~, ~, V] = svd(H);

for i = 1 : Ns
    V1(:,i) = V(:,i);
end

end