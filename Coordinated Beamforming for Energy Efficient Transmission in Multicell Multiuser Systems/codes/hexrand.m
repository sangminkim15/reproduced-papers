function [y] = hexrand(L, L0, BS, N)

% I/O
% L - length of hexagon radius
% L0 - minimum distance between user and base station
% BS - location of base station
% N - number of users per base station
% y - location of users

y = zeros(N, 2);
i = 1;

while true
    x = [L * (2 * rand(1,1) - 1), sqrt(3) / 2 * L * (2 * rand(1,1) - 1)] + BS;
    
    if norm(x - BS) > L0
        if sqrt(3) * abs(x(1) - BS(1)) + abs(x(2) - BS(2)) < L * sqrt(3)
            y(i,:) = x;
            i = i + 1;
        end
    end
    
    if i > N
        break
    end

end

end