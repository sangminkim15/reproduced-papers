function [y] = G(h, W, alpha, variance, eta, xi)

% I/O
% W{j,1} - beamforming cell array [M(j) N(j)]
% alpha{j,1} - priority of user cell array, size [1 N(j)]
% variance{j,1} - circ symm complex Gaussian RN cell array, size [1 N(j)]
% eta - energyeff coefficient
% xi - inefficiency of power amplifier
% y - G (equivalent form)

% f1
y1 = f1(h, W, alpha, variance);

% norm(W{j,1}(:,k))
y2 = zeros(1,1);
K = size(W);
M = zeros(1,K(1));
N = zeros(1,K(1));

for i = 1 : K(1)
    [M(i),N(i)] = size(W{i,1});
    y2 = y2 + eta * xi * trace(W{i,1}'*W{i,1});
end

% G
y = y1 - y2;

end