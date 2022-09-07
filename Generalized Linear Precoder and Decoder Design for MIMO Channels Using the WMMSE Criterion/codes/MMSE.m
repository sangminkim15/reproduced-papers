function [] = MMSE ()

H22 = channel(2,2);
H42 = channel(2,4);
H62 = channel(2,6);

sigmadB = -3.9 : 0.1 : 0;
sigma = 10.^(sigmadB / 10);
SNR = (2 * sigma).^-1;
SNRdB = 10 * log(SNR) / log(10);

BER22 = zeros(size(sigmadB));
BER42 = zeros(size(sigmadB));
BER62 = zeros(size(sigmadB));

for i = 1 : length(sigmadB)
    [~, BER22(i)] = BER_MMSE(sigma(i), H22, 2);
    [~, BER42(i)] = BER_MMSE(sigma(i), H42, 2);
    [~, BER62(i)] = BER_MMSE(sigma(i), H62, 2);
end

semilogy(SNRdB, BER22, SNRdB, BER42, SNRdB, BER62, 'LineWidth', 2);
legend('Mt = 2', 'Mt = 4', 'Mt = 6');
xlabel('Total Transmit Power / Received Noise (dB)');
ylabel('Average BER');
xlim([-3 1]);
grid on;

end