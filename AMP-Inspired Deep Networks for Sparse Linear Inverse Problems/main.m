function [] = main()

% Scenario I - A drawn i.i.d. Gaussian

% Generate A, x, v, y = Ax + v
M = 250;
N = 500;
A = normrnd(0, M.^(-0.5), [M, N]);

p = 0.1;
x = binornd(1, p * ones(N, 1));
for i = 1 : N
    x(i) = x(i) * normrnd(0, 1);
end

sigma = norm(A * x);
v = normrnd(0, 0.01 * sigma , [M, 1]);

y = A * x + v;

% ISTA / FISTA
Xi = ista(A, y, 1e4, 5e-3, 0.1);
Xf = fista(A, y, 1e4, 5e-3, 0.1);

% NMSE calculation
iteration = zeros(1, 1e4);
NMSEi = zeros(1, 1e4);
NMSEf = zeros(1, 1e4);

for i = 1 : 1e4
    iteration(i) = i;
    NMSEi(i) = (norm(Xi(:,i) - x).^2) / (norm(x).^2);
    NMSEi(i) = 10 * (log(NMSEi(i)) / log(10));
    NMSEf(i) = (norm(Xf(:,i) - x).^2) / (norm(x).^2);
    NMSEf(i) = 10 * (log(NMSEf(i)) / log(10));
end

semilogx(iteration, NMSEi, iteration, NMSEf, 'LineWidth', 2);
legend('ISTA', 'FISTA');
xlabel('Iterations');
ylabel('NMSE (dB)');
grid on

end