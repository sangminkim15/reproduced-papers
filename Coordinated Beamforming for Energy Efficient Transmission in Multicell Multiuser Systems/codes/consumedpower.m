function [P] = consumedpower(Aj, hjj, alphaj, muj, sj, lambdaj)

% I/O
% Aj - jth row of cell array A{j,1}, A = Amatrix(h, alpha, mu, s)
% hjj - h{j,j}, size [M(j) N(j)]
% alphaj - jth row of cell array alpha{j,1}, size [1 N(j)]
% muj - jth row of cell array mu{j,1}, size [1 N(j)]
% sj - jth row of cell array s{j,1}, size [1 N(j)]
% lambdaj - Lagrange multiplier

% U, D, V - eigenvalue decomposition of A

% size def
K = size(Aj);   % [M(j),M(j)]
L = size(hjj);  % [M(j),N(j)]

% eigenvalue decomposition
[U, D, V] = eig(Aj);

% PSI
PSI = zeros(K);

for k = 1 : L(2)
    PSI = PSI + norm(alphaj(k) * sj(k) * muj(k))^2 * hjj(:,k) * hjj(:,k)';
end
PSI = V' * PSI * U;

% power constraint
P = 0;

for m = 1 : K(1)
    P = P + (PSI(m,m) / (D(m,m) + lambdaj)^2);
end

end