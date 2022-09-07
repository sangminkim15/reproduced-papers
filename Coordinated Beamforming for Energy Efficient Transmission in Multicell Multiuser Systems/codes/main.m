function [] = main()

% CHANNEL

L = 0.5;
L0 = 0.4;
BS = cell(1,3);
BS{1} = [-0.5 0];
BS{2} = [0.25 0.25*sqrt(3)];
BS{3} = [0.25 -0.25*sqrt(3)];

d = distance(L, L0, BS, 3, [3 3 3]);
h = channel(d, [4 4 4]);

% PLOT

PdBm = 26 : 4 : 46;
P = 10.^(PdBm/10)/1000;
y = zeros(1, length(P));

for i = 1 : length(P)
    y(i) = energyeff(P(i), h);
end

plot(PdBm, y, '-*' , 'LineWidth', 2);
xlabel('Power Constraint (dB)');
ylabel('Energy Efficiency');
grid on;

end