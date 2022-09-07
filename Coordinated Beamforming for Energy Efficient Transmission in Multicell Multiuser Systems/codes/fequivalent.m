function [y] = fequivalent(h, W, alpha, variance, eta, Pc, P0, xi)

y = f1(h, W, alpha, variance) - eta * f2(W, Pc, P0, xi);

end