function [y] = f2(W, Pc, P0, xi)

% I/O
% W{j,1} - beamforming cell array [W(j,:,1) W(j,:,2) ... W(j,:,N(j))]
% Pc - constant circuit power per antenna
% P0 - basic power consumed at BS
% y - f2 (total power consumption)

% M(i) - # of antennas in cell i
% N(i) - # of users in cell i

y = zeros(1,1);
K = size(W);
M = zeros(1,K(1));
N = zeros(1,K(1));

for i = 1 : K(1)
    [M(i),N(i)] = size(W{i,1});
    y = y + xi * trace(W{i,1}'*W{i,1}) + M(i)*Pc + P0;
end

end