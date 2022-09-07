function [W2] = innerlayer(h, alpha, variance, eta, xi, Px, delta)

K = size(h);

W1 = cell(K(1), 1);
W2 = cell(K(1), 1);
mu = cell(K(1), 1);
e = cell(K(1), 1);
s = cell(K(1), 1);

L = zeros(1, K(1));
N = zeros(1, K(1));

for j = 1 : K(1)
    [L(j),N(j)] = size(h{j,j});
end

% STEP 1 : Initialize W{j,1}, mu{j,1}(k), s{j,1}(k)
n = 0;
for j = 1 : K(1)
    X = size(h{j,j});
    W1{j,1} = ones(size(h{j,j})) * (sqrt(Px(j)) / sqrt(X(2) * X(1)) / 10);
    W2{j,1} = ones(size(h{j,j})) * (sqrt(Px(j)) / sqrt(X(2) * X(1)) / 10);
    
    mu{j,1} = zeros(1, N(j));
    e{j,1} = zeros(1, N(j));
    s{j,1} = zeros(1, N(j));
end

% mu, e, s initialized to zero

% STEP 2 ~ 5
while true
    % STEP 2 : n = n + 1
    n = n + 1;
    
    W1 = W2;
    
    % STEP 3 : update mu, e, s using HH
    [mu, ~, s] = HH(h, W1, variance);
    
    % STEP 4 : update W2{j,1}(:,k) - use Amatrix -> consumedpower -> beamforming
    A = Amatrix(h, alpha, mu, s, eta, xi);
    
    for j = 1 : K(1)
        W2{j,1} = beamforming(A{j,1}, h{j,j}, alpha{j,1}, mu{j,1}, s{j,1}, Px(j));
    end
    
    % STEP 5. if abs(W2{j,1}(k) - W1{j,1}(k)) < delta -> break
    if abs(G(h, W2, alpha, variance, eta, xi) - G(h, W1, alpha, variance, eta, xi)) < delta
        break
    end

end

end