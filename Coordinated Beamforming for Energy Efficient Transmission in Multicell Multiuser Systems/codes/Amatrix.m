function [A] = Amatrix(h, alpha, mu, s, eta, xi)

% I/O
% h{m,j} - channel cell array, size [M(m) N(j)]
% alpha{j,1} - priority of user cell array, size [1 N(j)]
% mu{j,1} - optimal receiver cell array [1 N(j)]
% s{j,1} - auxiliary variables cell array [1 N(j)]
% A{j,1} - cell array, size [M(j) M(j)]

K = size(h);
M = zeros(1, K(1));
N = zeros(1, K(1));

% size def
for m = 1 : K(1)
    [M(m), N(m)] = size(h{m,m});
end

A = cell(K(1), 1);

% A{j,1}

for j = 1 : K(1)
    A{j,1} = zeros(M(j), M(j));
    for m = 1 : K(1)
        for n = 1 : N(m)
            A{j,1} = A{j,1} + alpha{m,1}(n) * s{m,1}(n) * norm(mu{m,1}(n))^2 * h{j,m}(:,n) * h{j,m}(:,n)';
        end
    end
    A{j,1} = A{j,1} + eta * xi * eye(M(j));
end

end