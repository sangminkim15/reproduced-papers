function [W] = outerlayer(h, alpha, variance, xi, Pc, P0, Px, etamin, etamax, epsilon, delta)

% Outer Layer Solution

while abs(etamax - etamin) > epsilon
    eta = (etamin + etamax) / 2;
    
    % Sub-Problem Solution { }
    % innerlayer.m 
    W = innerlayer(h, alpha, variance, eta, xi, Px, delta);
    
    y = fequivalent(h, W, alpha, variance, eta, Pc, P0, xi);
    if y > 0
        etamin = eta;
    else
        etamax = eta;
    end
end

end