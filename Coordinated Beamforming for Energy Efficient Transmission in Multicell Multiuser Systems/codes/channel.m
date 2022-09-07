function [h] = channel(d, M)

% I/O
% d{m,j} - distance (between user and BS) cell array, size [1, N(j)]
% M(m) - # of antennas in cell m
% h{m,j} - channel cell array, size [M(m), N(j)]

% hlarge{m,j} - large scaling fading
% hsmall{m,j} - small scaling fading

K = size(d);
L = zeros(1, K(1));
N = zeros(1, K(1));

% size def
for m = 1 : K(1)
    [L(m), N(m)] = size(d{m,m});
end

% initialize cell arrays
hlarge = cell(K(1), K(1));
hsmall = cell(K(1), K(1));
h = cell(K(1), K(1));

% channel
for j = 1 : K(1)
    for k = 1 : N(j)
        for m = 1 : K(1)
            hlarge{m,j}(k) = 10.^((-34.5 - 38 * log(d{m,j}(k))/log(10) + normrnd(0,8)) / 10);
            hsmall{m,j}(:,k) = normrnd(0,1/sqrt(2),[M(m),1]) + 1i * normrnd(0,1/sqrt(2),[M(m),1]);
            
            h{m,j}(:,k) = hlarge{m,j}(k) * eye(M(m)) * hsmall{m,j}(:,k);
        end
    end
end

end