clear
clc
close all

%% Intelligent Reflecting Surface-Aided Joint Processing Coordinated Multipoint Transmission
% IEEE Trans. Commun. available at https://ieeexplore.ieee.org/document/9279253/

%% Common parameters
p.cell_type = 'hexagonal';
p.side_length = 200*sqrt(3);
p.num_cell = 1;
p.simulation_scenario = 1;
p.L0_dB = -30; % channel power gain [dB]
p.L0 = 10^(p.L0_dB/10);
p.d0 = 1; %reference distance [m]
% p.dx = ; % link distance
p.alpha_br = 2.2; % path loss exponent of BS-IRS link (Rician)
p.alpha_ru = 2.2; % path loss exponent of IRS-user link (Rician)
p.alpha_bu = 3.6; % path loss exponent of BS-user link (Rayleigh)
p.Rician_f_dB = 10; % Rician factor [dB]
p.Rician_f = 10^(p.Rician_f_dB/10);
p.range_AoA = [0,2*pi]; % arrival of angle is randomly distributed within [0; 2*pi]
p.range_AoD = [0,2*pi]; % arrival of departure is randomly distributed within [0; 2*pi]
p.d = 2;   % number of desired data streams
p.N_r = 2; % number of reciver antennas at user
p.variance = -80; % noise variance [dB]
p.np = 10^((p.variance-30)/10);
p.channel_realization = 1; % channel realizations
p.b = 1; %  the number of bits to represent the resolution levels of IRS

p.num_BS = 2; % the number of BS
p.BS_1_location = [-300,0]; % [x coordinate y coordinate]
p.BS_2_location = [300,0];
p.cell_edge_user_1_location = [0,0];
p.IRS_location = [0,100*sqrt(3)];

p.epsilon = 1e-4;
p.zeta = 1e-6;
p.eta = 1e-5;

%% Specific parameters for Figure 2 Conversionce behaviour of Algorithm 3.
p.M = 100; % the number of IRS elements [20,50,100]
p.N_t = 2;

objective_value = [];
P_dB_range = 0 : 5 : 20;
temp = zeros(1,length(P_dB_range));

for exp_idx = 1:p.channel_realization
    H = channel_realization(p);
    objective_value = [];
    disp(['# experiments : ',num2str(exp_idx)])
    for p_idx = P_dB_range
        disp(['SNR : ',num2str(p_idx),'[dB]'])
        p.P_max_dB = p_idx;
        p.P_max = 10^(p.P_max_dB/10);
        objective_temp = real(algorithm_3(p,H));
        objective_value = [objective_value, objective_temp];
    end
    temp = temp + objective_value;
    
end
average_achievable_rate = temp/p.channel_realization;
figure
plot(P_dB_range,average_achievable_rate, 'LineWidth', 1.5)
grid on
title('Average achievable rate versus BS transmit power budget')
xlabel('Transmit power [dB]')
ylabel('Average achievable rate (bps/Hz)')
legend('MM, JP, continous')