function [mu, e, s] = HH(h, W, variance)

% I/O
% h{m,j} - channel cell array, size [M(m) N(j)]
% W{j,1} - beamforming cell array [M(j) N(j)]
% variance{j,1} - circ symm complex Gaussian RN cell array, size [1 N(j)]
% mu{j,1} - optimal receiver cell array [1 N(j)]
% s{j,1} - auxiliary variables cell array
% e{j,1} - MSE cell array [1 N(j)]

K = size(h);
M = zeros(1,K(1));
N = zeros(1,K(1));

% size def
for m = 1 : K(1)
    [M(m),N(m)] = size(h{m,m});
end

% initialize cell arrays
x = cell(K(1),1);
mu = cell(K(1),1);
e = cell(K(1),1);
s = cell(K(1),1);

for j = 1 : K(1)
    x{j,1} = zeros(1, N(j));
    mu{j,1} = zeros(1, N(j));
    e{j,1} = zeros(1, N(j));
    s{j,1} = zeros(1, N(j));
end

for j = 1 : K(1)
    for k = 1 : N(j)
        for m = 1 : K(1)
            x{j,1}(k) = x{j,1}(k) + h{m,j}(:,k)' * W{m,1} * W{m,1}' * h{m,j}(:,k);
        end
        mu{j,1}(k) = (h{j,j}(:,k)' * W{j,1}(:,k)) / (x{j,1}(k) + variance{j,1}(k));
        e{j,1}(k) = 1 - mu{j,1}(k);
        s{j,1}(k) = 1 / e{j,1}(k);
    end
end

end