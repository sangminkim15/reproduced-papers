function [] = EqualError ()

H = channel(5,5);

sigmadB = -3.9 : 0.1 : 0;
sigma = 10.^(sigmadB / 10);
SNR = (2 * sigma).^-1;
SNRdB = 10 * log(SNR) / log(10);

B3 = zeros(1, length(sigmadB));
B4 = zeros(1, length(sigmadB));
B5 = zeros(1, length(sigmadB));

for i = 1 : length(sigmadB)
    B3(i) = BER_EqualError(sigma(i) , H , 3);
    B4(i) = BER_EqualError(sigma(i) , H , 4);
    B5(i) = BER_EqualError(sigma(i) , H , 5);
end

semilogy(SNRdB, B3, SNRdB, B4, SNRdB, B5, 'LineWidth', 2);
legend('B = 3', 'B = 4', 'B = 5');
xlabel('Total Transmit Power / Received Noise (dB)');
ylabel('Average BER');
xlim([-3 1]);
grid on;

end