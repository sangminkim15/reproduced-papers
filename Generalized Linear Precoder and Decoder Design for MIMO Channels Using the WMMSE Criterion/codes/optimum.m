function [I, LAMBDA1, PHIf, PHIg] = optimum (B, LAMBDA, W, p0)

% Optimal Precoder - PHIf %

% I/O
% LAMBDA    Eigenvalue Matrix of H' * Rnn^(-1) * H
% W         Error Weight Matrix
% p0        Power Constraint

% mu        Lagrangian Multiplier
% PHIf      Optimum Precoder
% PHIg      Optimum Decoder

PHIf = zeros(B,B);
PHIg = zeros(B,B);

% Sort lambda * w in Descending Order
LAMBDA1 = zeros(size(LAMBDA));
W1 = zeros(size(W));

RHO = LAMBDA * W;

[RHO1, I] = sort(diag(RHO), 'descend');

for p = 1 : size(I)
    LAMBDA1(p,p) = LAMBDA(I(p),I(p));
    W1(p,p) = W(I(p),I(p));
end

% Lagrange Multiplier mu
k = B;          % Initialize k

while true
    x = 0;
    y = 0;
    
    for i = 1 : k
        x = x + LAMBDA1(i,i)^(-0.5) * W1(i,i)^0.5;
        y = y + LAMBDA1(i,i)^(-1);
    end
    
    mu = (x / (p0 + y))^2;
    
    if mu <= RHO1(k)
        for j = 1 : k
            PHIf(j,j) = mu^(-0.5) * LAMBDA1(j,j)^(-0.5) * W1(j,j)^0.5 - LAMBDA1(j,j)^(-1);
            
            if PHIf(j,j) >= 0
                PHIf(j,j) = PHIf(j,j)^0.5;
            else
                PHIf(j,j) = 0;
            end
        end
        
        break;
    else
        PHIf(k,k) = 0;
        k = k - 1;
    end
end

for l = 1 : B
    PHIg(l,l) = mu^0.5 * LAMBDA1(l,l)^(-0.5) * W1(l,l)^(-0.5) - mu * LAMBDA1(l,l)^(-1) * W1(l,l)^(-1);
    
    if PHIg(l,l) >= 0
        PHIg(l,l) = PHIg(l,l)^0.5 * LAMBDA1(l,l)^(-0.5);
    else
        PHIg(l,l) = 0;
    end
end
    
end