function [H] = channel(P, Q, K, SNR)

H = cell(K, 1);         % channel
for i = 1 : K
    H{i,1} = sqrt(SNR/2) * (randn(Q,P) + 1i * randn(Q,P));
end

end