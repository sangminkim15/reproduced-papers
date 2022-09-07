function [d] = distance(L, L0, BS, N, M)

% I/O
% L - length of hexagon diameter
% L0 - minimum distance between user and base station
% BS - location of base station, cell array
% N - number of base station
% M - number of users per base station
% d - distance between user and base station

% y - user location

% size def
y = cell(1, N);
d = cell(N, N);

% generation users
for i = 1 : N
    y{i} = hexrand(L, L0, BS{i}, M(i));
end

% distance
for m = 1 : N
    for j = 1 : N
        for k = 1 : M(j)
            d{m,j}(k) = norm(y{j}(k,:) - BS{m});
        end
    end
end

end