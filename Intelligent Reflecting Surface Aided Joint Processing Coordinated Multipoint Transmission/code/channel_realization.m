function [H] = channel_realization(p)

% Rician fading reference : https://arxiv.org/pdf/1907.10864.pdf

if p.simulation_scenario == 1 % Single Cell-Edge User
    
    % H.bs1_ue1 : H1
    % H.bs2_ue1 : H2
    % H.bs1_IRS : G1
    % H.bs2_IRS : G1
    % H.IRS_ue1 : Hr
    
    L_loss_bs1_ue1 = 10^(p.L0_dB/10)/(norm(p.BS_1_location-p.cell_edge_user_1_location)/p.d0)^(p.alpha_bu);
    L_loss_bs2_ue1 = 10^(p.L0_dB/10)/(norm(p.BS_2_location-p.cell_edge_user_1_location)/p.d0)^(p.alpha_bu);
    L_loss_bs1_IRS = 10^(p.L0_dB/10)/(norm(p.BS_1_location-p.IRS_location)/p.d0)^(p.alpha_br);
    L_loss_bs2_IRS = 10^(p.L0_dB/10)/(norm(p.BS_2_location-p.IRS_location)/p.d0)^(p.alpha_br); 
    L_loss_IRS_ue1 = 10^(p.L0_dB/10)/(norm(p.IRS_location-p.cell_edge_user_1_location)/p.d0)^(p.alpha_ru);
    
    % Rayleigh
    H.bs1_ue1 = (randn(p.N_r,p.N_t)+1i*randn(p.N_r,p.N_t))/sqrt(2)*sqrt(L_loss_bs1_ue1);
    % Rayleigh
    H.bs2_ue1 = (randn(p.N_r,p.N_t)+1i*randn(p.N_r,p.N_t))/sqrt(2)*sqrt(L_loss_bs2_ue1);
    % Rician
    H.bs1_IRS = Rician_fading(p,p.M,p.N_t)*sqrt(L_loss_bs1_IRS);
    % Rician
    H.bs2_IRS = Rician_fading(p,p.M,p.N_t)*sqrt(L_loss_bs2_IRS);
    % Rician
    H.IRS_ue1 = Rician_fading(p,p.N_r,p.M)*sqrt(L_loss_IRS_ue1);
    
else
    if p.simulation_scenario == 2 % Multiple Cell-Edge User
        
        L_loss_bs1_ue1 = 10^(p.L0_dB/10)/(norm(p.BS_1_location-p.cell_edge_user_1_location)/p.d0)^(p.alpha_bu);
        L_loss_bs2_ue1 = 10^(p.L0_dB/10)/(norm(p.BS_2_location-p.cell_edge_user_1_location)/p.d0)^(p.alpha_bu);
        L_loss_bs3_ue1 = 10^(p.L0_dB/10)/(norm(p.BS_3_location-p.cell_edge_user_1_location)/p.d0)^(p.alpha_bu);

        L_loss_bs1_ue2 = 10^(p.L0_dB/10)/(norm(p.BS_1_location-p.cell_edge_user_2_location)/p.d0)^(p.alpha_bu);
        L_loss_bs2_ue2 = 10^(p.L0_dB/10)/(norm(p.BS_2_location-p.cell_edge_user_2_location)/p.d0)^(p.alpha_bu);
        L_loss_bs3_ue2 = 10^(p.L0_dB/10)/(norm(p.BS_3_location-p.cell_edge_user_2_location)/p.d0)^(p.alpha_bu);
        
        L_loss_bs1_ue3 = 10^(p.L0_dB/10)/(norm(p.BS_1_location-p.cell_edge_user_3_location)/p.d0)^(p.alpha_bu);
        L_loss_bs2_ue3 = 10^(p.L0_dB/10)/(norm(p.BS_2_location-p.cell_edge_user_3_location)/p.d0)^(p.alpha_bu);
        L_loss_bs3_ue3 = 10^(p.L0_dB/10)/(norm(p.BS_3_location-p.cell_edge_user_3_location)/p.d0)^(p.alpha_bu);        
        
        L_loss_bs1_IRS = 10^(p.L0_dB/10)/(norm(p.BS_1_location-p.IRS_location)/p.d0)^(p.alpha_br);
        L_loss_bs2_IRS = 10^(p.L0_dB/10)/(norm(p.BS_2_location-p.IRS_location)/p.d0)^(p.alpha_br); 
        L_loss_bs3_IRS = 10^(p.L0_dB/10)/(norm(p.BS_3_location-p.IRS_location)/p.d0)^(p.alpha_br); 
        
        L_loss_IRS_ue1 = 10^(p.L0_dB/10)/(norm(p.IRS_location-p.cell_edge_user_1_location)/p.d0)^(p.alpha_ru);
        L_loss_IRS_ue2 = 10^(p.L0_dB/10)/(norm(p.IRS_location-p.cell_edge_user_2_location)/p.d0)^(p.alpha_ru);        
        L_loss_IRS_ue3 = 10^(p.L0_dB/10)/(norm(p.IRS_location-p.cell_edge_user_3_location)/p.d0)^(p.alpha_ru);
        
        % BS -> User 1
        H.bs1_ue1 = (randn(p.N_r,p.N_t)+1i*randn(p.N_r,p.N_t))/sqrt(2)*sqrt(L_loss_bs1_ue1);
        H.bs2_ue1 = (randn(p.N_r,p.N_t)+1i*randn(p.N_r,p.N_t))/sqrt(2)*sqrt(L_loss_bs2_ue1);
        H.bs3_ue1 = (randn(p.N_r,p.N_t)+1i*randn(p.N_r,p.N_t))/sqrt(2)*sqrt(L_loss_bs3_ue1);
        
        % BS -> User 2
        H.bs1_ue2 = (randn(p.N_r,p.N_t)+1i*randn(p.N_r,p.N_t))/sqrt(2)*sqrt(L_loss_bs1_ue2);
        H.bs2_ue2 = (randn(p.N_r,p.N_t)+1i*randn(p.N_r,p.N_t))/sqrt(2)*sqrt(L_loss_bs2_ue2);
        H.bs3_ue2 = (randn(p.N_r,p.N_t)+1i*randn(p.N_r,p.N_t))/sqrt(2)*sqrt(L_loss_bs3_ue2);
        
        % BS -> User 3
        H.bs1_ue3 = (randn(p.N_r,p.N_t)+1i*randn(p.N_r,p.N_t))/sqrt(2)*sqrt(L_loss_bs1_ue3);
        H.bs2_ue3 = (randn(p.N_r,p.N_t)+1i*randn(p.N_r,p.N_t))/sqrt(2)*sqrt(L_loss_bs2_ue3);
        H.bs3_ue3 = (randn(p.N_r,p.N_t)+1i*randn(p.N_r,p.N_t))/sqrt(2)*sqrt(L_loss_bs3_ue3);        
        
        % BS -> IRS
        H.bs1_IRS = Rician_fading(p,p.M,p.N_t)*sqrt(L_loss_bs1_IRS);
        H.bs2_IRS = Rician_fading(p,p.M,p.N_t)*sqrt(L_loss_bs2_IRS);
        H.bs3_IRS = Rician_fading(p,p.M,p.N_t)*sqrt(L_loss_bs3_IRS);
        
        % IRS -> Users
        H.IRS_ue1 = Rician_fading(p,p.N_r,p.M)*sqrt(L_loss_IRS_ue1);
        H.IRS_ue2 = Rician_fading(p,p.N_r,p.M)*sqrt(L_loss_IRS_ue2);
        H.IRS_ue3 = Rician_fading(p,p.N_r,p.M)*sqrt(L_loss_IRS_ue3);
        
    end
end

end