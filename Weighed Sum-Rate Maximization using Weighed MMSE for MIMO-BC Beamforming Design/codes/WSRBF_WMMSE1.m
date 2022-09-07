function [WSR] = WSRBF_WMMSE1 (u, H, Etx)

% Algorithm 1
%   10 random filter initializations &
%       select one that leads to the highest WSR

% I/O
% u     Weight, size 1 * K
% H     Channel Matrix, size QK * P
% Etx   Power Constraint

% WSR   Weight Sum Rate

% Initialization
[K, ~] = size(H);
[Q, P] = size(H{1,1});

B = cell(1,K);

A = cell(K,K);
W = cell(K,K);

E = cell(1,K);
R = zeros(1,K);

x = zeros(1,10);
y = zeros(1,10);

% 10 Random Initialization
n = 0;
while n < 10
    tr = 0;
    for i = 1 : K
        B{1,i} = sqrt(1/2) * (randn(P,Q) + 1i * randn(P,Q));
        tr = tr + trace(B{1,i} * B{1,i}');
    end
    
    for i = 1 : K
        B{1,i} = (Etx / tr).^0.5 * B{1,i};
    end
    
    while true
        y(1,n+1) = x(1,n+1);
        
        for j = 1 : K
            for l = 1 : K
                if j == l
                    A{j,l} = RxFilter(H{j,1}, B, j);
                    [~, W{j,l}] = WMMSE(u(j), H{j,1}, B, j);
                else
                    A{j,l} = zeros(Q,Q);
                    W{j,l} = zeros(Q,Q);
                end
            end
        end
    
        Hmatrix = cell2mat(H);
        Amatrix = cell2mat(A);
        Wmatrix = cell2mat(W);
    
        Bmatrix = Beamforming(Hmatrix, Amatrix, Wmatrix, Etx);
    
        B = matrixdiv(Bmatrix, K);
    
        for p = 1 : K
            [E{1,p}, ~] = WMMSE(u(p), H{p,1}, B, p);
            R(p) = log(det(E{1,p}^(-1))) / log(2);
        end
        
        x(1,n+1) = sum(u.*R);
        
        if abs(x(1,n+1) - y(1,n+1)) < 1e-2
            n = n + 1;
            break
        end
    end
end

WSR = max(x);

end