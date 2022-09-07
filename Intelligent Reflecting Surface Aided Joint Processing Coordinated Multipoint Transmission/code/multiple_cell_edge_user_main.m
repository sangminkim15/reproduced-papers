clear all,clc,close all

%% Intelligent Reflecting Surface-Aided Joint Processing Coordinated Multipoint Transmission
% IEEE Trans. Commun. available at https://ieeexplore.ieee.org/document/9279253/

%% Common parameters
p.cell_type = 'hexagonal';
p.side_length = 200*sqrt(3);
p.num_cell = 1;
p.simulation_scenario = 2; % 1-single cell-edge user , 2-Multiple cell-edge user
p.L0_dB = -30; % channel power gain [dB]
p.L0 = 10^(p.L0_dB/10);
p.d0 = 1; %reference distance [m]
p.num_BS = 3; % the number of BS
p.BS_1_location = [-300,0]; % [x coordinate y coordinate]
p.BS_2_location = [300,0];
p.BS_3_location = [0,300*sqrt(3)];
p.IRS_location = [0,100*sqrt(3)];
p.central_point_cell_edge_users = [0,100*sqrt(3)];
p.circle_radius = 30; % [m]
p.IRS_altitude = 10; % [m]
p.cell_edge_user_1_location = p.IRS_location + p.circle_radius .* [cos(2*pi*rand(1,1)), sin(2*pi*rand(1,1))];
p.cell_edge_user_2_location = p.IRS_location + p.circle_radius .* [cos(2*pi*rand(1,1)), sin(2*pi*rand(1,1))];
p.cell_edge_user_3_location = p.IRS_location + p.circle_radius .* [cos(2*pi*rand(1,1)), sin(2*pi*rand(1,1))];
p.N_t = 6; % the number of transmit antennas
p.N_r = 6; % the number of received antennas
p.M = 100; % the number of IRS elements [20,50,100]
p.b = 2; %  the number of bits to represent the resolution levels of IRS
p.K = 3; % the number of users

p.alpha_br = 2.2; % path loss exponent of BS-IRS link (Rician)
p.alpha_ru = 2.2; % path loss exponent of IRS-user link (Rician)
p.alpha_bu = 3.6; % path loss exponent of BS-user link (Rayleigh)
p.Rician_f_dB = 10; % Rician factor [dB]
p.Rician_f = 10^(p.Rician_f_dB/10);
p.range_AoA = [0,2*pi]; % arrival of angle is randomly distributed within [0; 2*pi]
p.range_AoD = [0,2*pi]; % arrival of departure is randomly distributed within [0; 2*pi]
p.d = 2;   % number of desired data streams
p.variance = -80; % noise variance [dB]
p.np = 10^(p.variance/10) * 0.001;
p.channel_realization = 10; % channel realizations

p.eta = 1e-3;


%% Fig.9 Average max-min rate versus BS tranmit power budget
objective_value = [];
P_dB_range = -10 : 5 : 10;
temp = zeros(1,length(P_dB_range));

for exp_idx = 1:p.channel_realization
    H = channel_realization(p);
    objective_value = [];
    disp(['# experiments : ',num2str(exp_idx)])
    for p_idx = P_dB_range
        disp(['SNR : ',num2str(p_idx),'[dB]'])
        p.P_max_dB = p_idx;
        p.P_max = 10^(p.P_max_dB/10);
        % Multiple cell-edge users (MCEU)
        objective_temp = algorithm_MCEU(p,H); 
        objective_value = [objective_value, objective_temp];
    end
    temp = temp + objective_value;
    
end
avg_max_min_rate = temp/p.channel_realization;
figure
plot(P_dB_range,avg_max_min_rate, 'LineWidth', 1.5)
grid on
title('Average max-min rate versus BS transmit power budget')
xlabel('Transmit power [dB]')
ylabel('Average achievable rate [bps/Hz]')
legend('Optimized phase, JP, continous','Optimized phase, JP, 2-bit','Optimized phase, JP, 1-bit')