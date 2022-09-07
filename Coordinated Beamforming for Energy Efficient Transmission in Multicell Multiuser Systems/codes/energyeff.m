function [y] = energyeff(P, h)

% ENERGY EFFICIENCY

% PARAMETER SETTINGS
Pc = 1;     % 30dBm
P0 = 10;    % 40dBm
xi = 1;
alpha = cell(3,1);
for i = 1 : 3
    alpha{i,1} = [1, 1, 1];
end

% VARIANCE SETTINGS
variance = cell(3,1);
for i = 1 : 3
    variance{i,1} = [10.^(-0.9), 10.^(-0.9), 10.^(-0.9)];
end

% ERROR SETTINGS
epsilon = 1e-5;
delta = 1e-3;

% POWER CONSTRAINT SETTINGS
Px = ones(1,3) * P;

% OPTIMAL BEAMFORMING
W = outerlayer(h, alpha, variance, xi, Pc, P0, Px, 0, 1e10, epsilon, delta);
y = f1(h, W, alpha, variance) / f2(W, Pc, P0, xi);

end