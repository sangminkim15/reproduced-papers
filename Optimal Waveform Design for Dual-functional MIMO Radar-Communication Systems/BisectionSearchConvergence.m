function [X] = BisectionSearchConvergence (Q, G, lambda_low, lambda_high, p)

Norm = p.L * p.Pt;
lambda = (lambda_low + lambda_high) / 2;

X = (Q + lambda * eye(p.N, p.N))^(-1) * G;
NormTemp = norm(X, 'fro')^2;

idx = 0;
idxarray = 0;
Normarray = NormTemp;

while abs(NormTemp - Norm) >= 1e-6
    idx = idx + 1;
    if NormTemp < Norm
        lambda_high = lambda;
        lambda = (lambda_low + lambda_high) / 2;
        
        X = (Q + lambda * eye(p.N, p.N))^(-1) * G;
        NormTemp = norm(X, 'fro')^2;
    
    else
        lambda_low = lambda;
        lambda = (lambda_low + lambda_high) / 2;
        
        X = (Q + lambda * eye(p.N, p.N))^(-1) * G;
        NormTemp = norm(X, 'fro')^2;
        
    end
    
    idxarray = horzcat(idxarray, idx);
    Normarray = horzcat(Normarray, NormTemp);
    
end

plot(idxarray, Normarray, 'LineWidth', 1.5);
title('');
xlabel('# of iterations');
ylabel('|X(\lambda)|^2');

end