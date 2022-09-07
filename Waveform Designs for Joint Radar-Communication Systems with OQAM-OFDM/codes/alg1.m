function [S] = alg1 (SNRin, K, Rrad, rho, Q, d)

P = zeros(K, Rrad);
p = zeros(K, 1);
for idx = 1 : Rrad
    P(:, idx) = circshift(Q(:, idx), Rrad-idx);
    p = p + P(:, idx);
end
p = flipud(p);              
p = (1/K) * fft(p);         % p = [p(0) p(1) ... p(K-1)]
x = p / norm(p);            % x = [x(0) x(1) ... x(K-1)]

d_theta = zeros(K, 1);
% angle optimization
for jdx = 1 : K
    d_theta(jdx, 1) = exp(1i * angle(d(jdx, 1)));
end
    
% magnitude optimization
x_abs = abs(x);
x_norm = 1;
d_norm = SNRin;
remove_idx_array = zeros(0,0);
d_abs = x_abs .* sqrt(d_norm / x_norm);
while true
    logic = d_abs >= rho;
    if logic == ones(K, 1)
        break
    
    else
        %%%%% remove d_abs s.t. smaller than 0 and closest to 0. %%%%%
        r = zeros(0, 0);
        for kdx = 1 : K
            if logic(kdx) == 0
                r = horzcat(r, [d_abs(kdx)-rho(kdx) ; kdx]);
            end
        end
        
        remove_num = find(r(1,:) == max(r(1,:)));
        remove_idx = r(2,remove_num(1));
        
        remove_idx_array = horzcat(remove_idx_array, remove_idx);
        
        d_abs(remove_idx) = rho(remove_idx);
        x_norm = x_norm - x_abs(remove_idx).^2;
        d_norm = d_norm - rho(remove_idx).^2;
        
        for ldx = 1 : K
            if any(remove_idx_array == ldx)
                d_abs(ldx) = rho(ldx);
                
            else
                d_abs(ldx) = x_abs(ldx) .* sqrt(d_norm / x_norm);
            end
        end
        
    end

end

d = d_abs .* d_theta;

s = K * ifft(d);       
S = zeros(K, Rrad);
for mdx = 1 : Rrad
    S(:, mdx) = circshift(flipud(s), K-Rrad+mdx);
end

end