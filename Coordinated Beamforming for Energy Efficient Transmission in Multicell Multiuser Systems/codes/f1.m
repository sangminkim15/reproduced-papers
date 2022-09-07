function [y] = f1(h, W, alpha, variance)

% I/O
% h{m,j} - channel cell array, size [M(m) N(j)]
% W{j,1} - beamforming cell array [M(j) N(j)]
% alpha{j,1} - priority of user cell array, size [1 N(j)]
% variance{j,1} - circ symm complex Gaussian RN cell array, size [1 N(j)]
% y - f1 (weighted sum rate)

% M(m) - # of antennas in cell m
% N(m) - # of users in cell m
% gamma{j,1} - interference cell array, size [1 N(j)]
% R{j,1} - rate cell array, size [1 N(j)]

K = size(W);
M = zeros(1,K(1));
N = zeros(1,K(1));

% size def
for m = 1 : K(1)
    [M(m),N(m)] = size(W{m,1});
end

% initialize cell arrays
x = cell(K(1),1);
gamma = cell(K(1),1);
R = cell(K(1),1);

for j = 1 : K(1)
    x{j,1} = zeros(1,N(j));
    gamma{j,1} = zeros(1,N(j));
    R{j,1} = zeros(1,N(j));
end

y = zeros(1,1);

% f1
for j = 1 : K(1)
    for k = 1 : N(j)
        for m = 1 : K(1)
            x{j,1}(k) = x{j,1}(k) + h{m,j}(:,k)' * W{m,1} * W{m,1}' * h{m,j}(:,k);
        end
        gamma{j,1}(k) = x{j,1}(k) - norm(h{j,j}(:,k)' * W{j,1}(:,k))^2;
        R{j,1}(k) = log(1 + norm((h{j,j}(:,k)' * W{j,1}(:,k)))^2/(gamma{j,1}(:,k) + variance{j,1}(k)));
        y = y + alpha{j,1}(k) * R{j,1}(k);
    end
end

end