function [y] = ArrayResponse_vec(N, phi, theta)

% Generates Normalized Tx and Rx Array Vectors

% I/O
% N         square root of # of transmit antennas (Nt = N^2)
% phi       azimuth angle
% theta     elevation angle
% y         normalized transmit array received vector

y = zeros(N^2, 1);

for i = 1 : N
    for j = 1 : N
        y((i-1) * N + j) = (1/N) * exp(1i * pi * ((i-1) * sin(phi) * sin(theta) + (j-1) * cos(theta)));
    end
end

end