function [H] = Rician_fading(p,N_r,N_t)

theta_AoA = exp(j*2*pi*(1/2)*[0:(N_r-1)]'*sin(rand(1)*2*pi));
theta_AoD = exp(j*2*pi*(1/2)*[0:(N_t-1)]'*sin(rand(1)*2*pi));

H_Los = theta_AoA*theta_AoD';
H_NLoS = (randn(N_r,N_t)+j*randn(N_r,N_t))/sqrt(2);

H = sqrt(p.Rician_f/(p.Rician_f+1))*H_Los + sqrt(1/(p.Rician_f+1))*H_NLoS;
end