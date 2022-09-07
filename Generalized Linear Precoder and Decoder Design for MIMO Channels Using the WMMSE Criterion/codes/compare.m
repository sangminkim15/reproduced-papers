function [] = compare ()

% Comparison of MMSE Design and Equal Error Design % 

H = channel(4,4);

sigmadB = -3.9 : 0.1 : 0;
sigma = 10.^(sigmadB / 10);
SNR = (2 * sigma).^-1;
SNRdB = 10 * log(SNR) / log(10);

MMSE1 = zeros(1,length(sigmadB));
MMSE2 = zeros(1,length(sigmadB));
MMSE3 = zeros(1,length(sigmadB));
EqualError = zeros(1,length(sigmadB));

for i = 1 : length(sigmadB)
    [MMSE, ~] = BER_MMSE(sigma(i), H, 3);
    MMSE1(i) = MMSE(1,1);
    MMSE2(i) = MMSE(2,2);
    MMSE3(i) = MMSE(3,3);
    
    EqualError(i) = BER_EqualError(sigma(i), H, 3);
end

semilogy(SNRdB, MMSE1, SNRdB, MMSE2, SNRdB, MMSE3, SNRdB, EqualError, 'LineWidth', 2);
legend('MMSE 1', 'MMSE 2', 'MMSE 3', 'Equal Error');
xlabel('Total Transmit Power / Received Noise (dB)');
ylabel('Average BER');
xlim([-3 1]);
grid on;

end