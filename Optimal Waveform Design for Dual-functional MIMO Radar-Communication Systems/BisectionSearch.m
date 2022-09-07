function [X] = BisectionSearch (Q, G, lambda_low, lambda_high, p)

Norm = p.L * p.Pt;
lambda = (lambda_low + lambda_high) / 2;

X = (Q + lambda * eye(p.N, p.N)) \ G;
NormTemp = norm(X, 'fro')^2;

while abs(NormTemp - Norm) >= 1e-6
    if NormTemp < Norm
        lambda_high = lambda;
        lambda = (lambda_low + lambda_high) / 2;
        
        X = (Q + lambda * eye(p.N, p.N)) \ G;
        NormTemp = norm(X, 'fro')^2;
    
    else
        lambda_low = lambda;
        lambda = (lambda_low + lambda_high) / 2;
        
        X = (Q + lambda * eye(p.N, p.N)) \ G;
        NormTemp = norm(X, 'fro')^2;
        
    end
    
end

end