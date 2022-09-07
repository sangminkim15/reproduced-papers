function [] = QoS ()

H = channel(3,3);
D = 1 / (10^0.5 + 1) * [10^0.5 0 ; 0 1];

sigmadB = -3.9 : 0.1 : 0;
sigma = 10.^(sigmadB / 10);
SNR = (2 * sigma).^-1;
SNRdB = 10 * log(SNR) / log(10);

SNR1 = zeros(1, length(sigmadB));
SNR2 = zeros(1, length(sigmadB));

for i = 1 : length(sigmadB)
    [SNR1(i), SNR2(i)] = SNR_QoS(sigma(i), H, D);
end

semilogy(SNRdB, SNR1, SNRdB, SNR2, 'LineWidth', 2)
legend('Video Stream', 'Audio Stream');
xlabel('Total Transmit Power / Received Noise (dB)');
ylabel('Sub Channel SNR');
xlim([-3 1]);
grid on;

end