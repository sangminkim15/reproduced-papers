function [Wj] = beamforming(Aj, hjj, alphaj, muj, sj, Px)

% I/O
% hjj - h{j,j}
% alphaj - jth row of cell array alpha{j,1}, size [1 N(j)]
% muj - jth row of cell array mu{j,1}, size [1 N(j)]
% sj - jth row of cell array s{j,1}, size [1 N(j)]
% Px - power constraint of Wj (scalar)
% Wj - jth row of beamforming cell array W{j,1}, size[M(j) N(j)]

lambdamin = 0;
lambdamax = 1e10;
lambdaj = 0;
epsilon = 1e-5;

x = consumedpower(Aj, hjj, alphaj, muj, sj, 0);

while abs(lambdamax - lambdamin) > epsilon
    if x <= Px
        break;
    else
        lambdaj = (lambdamin + lambdamax) / 2;
        y = consumedpower(Aj, hjj, alphaj, muj, sj, lambdaj);
        if y > Px
            lambdamin = lambdaj;
        else
            lambdamax = lambdaj;
        end
    end
end

% size def
K = size(Aj);   % [M(j),M(j)]
L = size(hjj);  % [M(j),N(j)]
Wj = zeros(L);

for k = 1 : L(2)
    Wj(:,k) = alphaj(k) * sj(k) * muj(k) * pinv(Aj + lambdaj * eye(K(1))) * hjj(:,k);
end

end