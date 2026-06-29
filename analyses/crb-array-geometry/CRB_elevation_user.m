clear; close all; clc;
tic;

%% Parâmetros Gerais e do Cenário
dd      = 50;
scale   = 1.5;
L       = 10;
x_min   = -50; x_max = 50;
y_min   = 10;  y_max = 50;
Nsub    = 1;
Nx = 8; Ny = 8;
dspace     = 0.5;
f_c        = 15e9;
lambda     = 3e8 / f_c;
d          = 0.5 * lambda;
P_tx       = 1;
SNR_dB     = 10;
SNR        = 10^(SNR_dB/10);

%% Definição das posições do usuário
user_positions = [30, 30, 0; 
                  0, 30*sqrt(2), 0;
                  15*sqrt(7), 15, 0;];
num_users = size(user_positions,1);

%% Vetor de alturas (0 a 100, 1000 pontos)
z_vec = linspace(0,100,1000);
CRB_norm_all = zeros(num_users,length(z_vec));

%% Loop sobre as posições do usuário
for u = 1:num_users
    UEs = user_positions(u,:);
    x_ue = UEs(1);
    y_ue = UEs(2);
    
    for i = 1:length(z_vec)
        z_antenna = z_vec(i);
        
        % Gera posições do array 8x8 com altura z_antenna
        elements_x_a = zeros(Nsub*Nx, Ny);
        elements_y_a = zeros(Nsub*Nx, Ny);
        elements_z_a = zeros(Nsub*Nx, Ny);
        x_offset = ((Nx-1)*d/2) + (Nsub-1)*(Nx*d + dspace)/2;
        z_offset = ((Ny-1)*d/2);
        for k_sub = 1:Nsub
            for i_ant = 1:Nx
                for j_ant = 1:Ny
                    idx = (k_sub-1)*Nx + i_ant;
                    elements_x_a(idx, j_ant) = (i_ant-1)*d + (k_sub-1)*(Nx*d + dspace) - x_offset;
                    elements_y_a(idx, j_ant) = 0;
                    elements_z_a(idx, j_ant) = (j_ant-1)*d - z_offset + z_antenna;
                end
            end
        end
        elements_x = elements_x_a(:);
        elements_y = elements_y_a(:);
        elements_z = elements_z_a(:);
        
        % Referência
        x_k1 = elements_x(1);
        y_k1 = elements_y(1);
        z_k1 = elements_z(1);
        
        % Cálculo dos termos para a FIM
        delta_1 = abs(x_ue - elements_x(1));
        r_1 = sqrt((x_ue-elements_x(1))^2 + y_ue^2 + elements_z(1)^2);
        delta_n = abs(x_ue - elements_x);
        r_n = sqrt((x_ue-elements_x).^2 + y_ue^2 + elements_z.^2);
        
        Term1_Jxx = (delta_n.^2) ./ (r_n.^6);
        Term2_Jxx = (2*pi/lambda)^2 * ( delta_n./(r_n.^2) - delta_1./(r_1*r_n) ).^2;
        Jxx = sum(Term1_Jxx + Term2_Jxx);
        
        Term1_Jyy = y_ue^2 ./ (r_n.^6);
        Term2_Jyy = (2*pi/lambda)^2 * ( 1./(r_n.^2) - 1./(r_1*r_n) ).^2;
        Jyy = sum(Term1_Jyy + y_ue^2 * Term2_Jyy);
        
        Term1_Jxy = delta_n ./ (r_n.^6);
        Term2_Jxy = (2*pi/lambda)^2 * ((delta_n./(r_n.^2) - delta_1./(r_1*r_n)) .* (1./(r_n.^2) - 1./(r_1*r_n)));
        Jxy = sum( y_ue*(Term1_Jxy+Term2_Jxy) );
        
        % Cálculo da potência total recebida
        total_rx_power = 0;
        for n = 1:length(elements_x)
            d_n = sqrt((x_ue-elements_x(n))^2 + y_ue^2 + elements_z(n)^2);
            beta_n = (lambda/(4*pi*d_n))^2;
            total_rx_power = total_rx_power + P_tx*beta_n;
        end
        noise_power = total_rx_power / SNR;
        
        % Normalização e matriz de Fisher
        Pnorm = L * (lambda/(sqrt(8)*pi*sqrt(noise_power)))^2;
        Jxx_scaled = Pnorm * Jxx;
        Jyy_scaled = Pnorm * Jyy;
        Jxy_scaled = Pnorm * Jxy;
        
        % CRB em x e y e norma combinada
        CRB_den_x = Jxx_scaled - (Jxy_scaled^2 / Jyy_scaled);
        CRB_den_y = Jyy_scaled - (Jxy_scaled^2 / Jxx_scaled);
        CRB_x = 1/CRB_den_x;
        CRB_y = 1/CRB_den_y;
        CRB_norm = sqrt(CRB_x^2 + CRB_y^2);
        
        CRB_norm_all(u,i) = CRB_norm;
    end
end

%% Plot das curvas para cada posição do usuário
colors = lines(num_users);
figure; hold on;
for u = 1:num_users
    plot(z_vec, CRB_norm_all(u,:), 'LineWidth',2, 'Color', colors(u,:), ...
         'DisplayName', sprintf('User (%.2f, %.2f, 0)', user_positions(u,1), user_positions(u,2)) );
end
xlabel('Antenna Mean Height (m)');
ylabel('Combined CRB');
title('Combined CRB vs. Antenna Height for Different User Positions');
grid on; legend('Location','best');
set(gca, 'YScale', 'log');

elapsedTime = toc;
fprintf('Elapsed time: %.2f s\n', elapsedTime);
