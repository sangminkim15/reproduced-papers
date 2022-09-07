function [H] = channel(Mr, Mt)

H = sqrt(1/2) * (randn(Mr,Mt) + 1i * randn(Mr,Mt));

end