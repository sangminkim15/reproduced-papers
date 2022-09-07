function [iterations, Zplot, SNRplot, SNRmaxplot, SNReqplot, ERRORplot, ERRORplot2, SS] = cyclic ()

% Parameters
K = 128;                        % # of subcarriers
Rcom = 6;                       % channel length
Rrad = 64;                      % # of range cells
SNRin = 10;                     % input SNR = 10dB
sigma = sqrt(1)/sqrt(K);        % std.

% Communication Channel
hdB = zeros(K, 1);
hdB(1:Rcom, 1) = [-6.0 0.0 -7.0 -22.0 -16.0 -20.0];    % power profile
h = zeros(K, 1);
h(1:Rcom, 1) = 10.^(hdB(1:Rcom, 1)/10);
H = fft(h);                                            % frequency response

% Threshold
SNRmin = 10.^(-20/10);                      % threshold SNR : -20dB per each subcarrier
rho = sigma * sqrt(SNRmin) ./ abs(H);       % threshold power of each subcarrier rho = [rho(0) rho(1) ... rho(K-1)]

% Initialization
d = randn(K, 1) .* ((2 .* randi([0 1], K, 1) - 1) + 1i * (2 .* randi([0 1], K, 1) - 1));
d = sqrt(SNRin) * d / norm(d);                % OFDM freq domain d = [d(0) d(1) ... d(K-1)]

dd = zeros(K, 1);
for idx = 1 : K
    dd(idx) = (2 .* randi([0 1], 1, 1) - 1) + 1i .* (2 .* randi([0 1], 1, 1) - 1);
    dd(idx) = dd(idx) / norm(dd(idx));
end
dd = sqrt(SNRin) * dd / norm(dd);

ss = K * ifft(dd);
S_S = zeros(K, Rrad);
for idx = 1 : Rrad
    S_S(:,idx) = circshift(flipud(ss), K-Rrad+idx);
end

Rand = orth(randn(K, K));
Q = sqrt(K)* sqrt(SNRin) * Rand(:, 1:Rrad);            % semiunitary matrix
s = K * ifft(d);                                       % OFDM time domain s = [s(0) s(1) ... s(K-1)]
S = zeros(K, Rrad);
for idx = 1 : Rrad
    S(:,idx) = circshift(flipud(s), K-Rrad+idx);
end

% Objective Function
y = norm(S'*S - (Q'* Q), 'fro');

jdx = 0;
z = y;

iterations = zeros(0, 0);
Zplot = zeros(0, 0);
SNRplot = zeros(0, 0);
SNRmaxplot = zeros(0, 0);
SNReqplot = zeros(0, 0);
ERRORplot = zeros(0, 0);
ERRORplot2 = zeros(0, 0);

while true
    y = z;
    % #1 : Obtain S given Q
    S = alg1 (SNRin, K, Rrad, rho, Q, d);
    
    % #2 : Obtain Q given S
    [U, ~, V] = svd(S');
    Vtilde = V(:, 1:Rrad);
    Q = sqrt(K) * sqrt(SNRin) * Vtilde * U';
    
    z = norm(S'*S - (Q'*Q), 'fro');
    
    SS = S' * S;
    
    SNRrad = sum(1 ./ diag(inv(S'*S))) / Rrad;
    SNRraddB = 10 * log(SNRrad) / log(10);
    
    s = S(:,Rrad);
    d = (1/K) * fft(flipud(s));
    
    SNRmax = 1 / norm(s).^(-2);
    SNRmaxdB = 10 * log(SNRmax) / log(10);

    SNReq = sum(1 ./ diag(inv(S_S' * S_S))) / Rrad;
    SNReqdB = 10 * log(SNReq) / log(10);
    
    
    SNRcom = abs(d).^2 .* abs(H).^2 / sigma^2;
    ErrorP = sum(erfc(SNRcom ./ sqrt(2))) ./ K;
    
    SNRcom2 = abs(dd).^2 .* abs(H).^2 / sigma^2;
    ErrorP2 = sum(erfc(SNRcom2 ./ sqrt(2))) ./ K;
    
    jdx = jdx + 1;
    iterations = horzcat(iterations, jdx);
    Zplot = horzcat(Zplot, z);
    SNRplot = horzcat(SNRplot, SNRraddB);
    SNRmaxplot = horzcat(SNRmaxplot, SNRmaxdB);
    SNReqplot = horzcat(SNReqplot, SNReqdB);
    ERRORplot = horzcat(ERRORplot, ErrorP);
    ERRORplot2 = horzcat(ERRORplot2, ErrorP2);
    
    if abs(z-y) < 1e-5      % stopping criterion
        break
    end
    
end

end