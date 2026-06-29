clear; close all; clc;
tic;

%% Parâmetros Gerais
dd      = 50;
scale   = 1.5;
L       = 10; 
x_min   = -50; x_max = 50;
y_min   = 10;  y_max = 50;
x_lim   = dd * scale;
y_lim   = dd * scale;

%% Parâmetros Comuns do Cenário
Nsub       = 1;
dspace     = 0.5;
f_c        = 15e9;
lambda     = 3e8 / f_c;
d          = 0.5 * lambda;
P_tx       = 1;
UEs        = [30, 30, 0];  % posição do usuário
SNR_dB     = 10;           % SNR fixo
SNR        = 10^(SNR_dB/10);

%% Vetor de alturas (z de 0 a 100, 1000 pontos)
z_vec = linspace(0,100,1000);

%% Vetor dos valores de Nx e Ny (Nx = Ny)
array_sizes = 2:10;
N_configs = length(array_sizes);

%% Matrizes para armazenar os resultados para cada configuração
CRB_norm_all = zeros(N_configs, length(z_vec));
CRB_x_all    = zeros(N_configs, length(z_vec));
CRB_y_all    = zeros(N_configs, length(z_vec));

%% Loop sobre diferentes configurações de antenas
for idx = 1:N_configs
    % Para cada configuração, Nx = Ny = current_size
    N_ant = array_sizes(idx);
    Nx = N_ant;
    Ny = N_ant;
    
    % Loop sobre altura das antenas
    for i = 1:length(z_vec)
        z_antenna = z_vec(i);
        
        % Gera posições das antenas com a altura z_antenna
        elements_x_a = zeros(Nsub*Nx, Ny);
        elements_y_a = zeros(Nsub*Nx, Ny);
        elements_z_a = zeros(Nsub*Nx, Ny);
        x_offset = ((Nx-1)*d/2) + (Nsub-1)*(Nx*d + dspace)/2;
        z_offset = ((Ny-1)*d/2);
        for k_sub = 1:Nsub
            for i_ant = 1:Nx
                for j_ant = 1:Ny
                    index = (k_sub-1)*Nx + i_ant;
                    elements_x_a(index, j_ant) = (i_ant-1)*d + (k_sub-1)*(Nx*d + dspace) - x_offset;
                    elements_y_a(index, j_ant) = 0;
                    elements_z_a(index, j_ant) = (j_ant-1)*d - z_offset + z_antenna;
                end
            end
        end
        
        % Vetores de posição dos elementos
        elements_x = elements_x_a(:);
        elements_y = elements_y_a(:);
        elements_z = elements_z_a(:);
        
        % Primeiro elemento como referência
        x_k1 = elements_x(1);
        y_k1 = elements_y(1);
        z_k1 = elements_z(1);
        
        % Cálculo dos termos para a FIM
        x_ue = UEs(1);
        y_ue = UEs(2);
        delta_1 = abs(x_ue - elements_x(1));
        r_1 = sqrt( (x_ue - elements_x(1))^2 + y_ue^2 + elements_z(1)^2 );
        delta_n = abs(x_ue - elements_x);
        r_n = sqrt( (x_ue - elements_x).^2 + y_ue^2 + elements_z.^2 );
        
        Term1_Jxx = (delta_n.^2) ./ (r_n.^6);
        Term2_Jxx = (2*pi/lambda)^2 * (delta_n./(r_n.^2) - delta_1./(r_1*r_n)).^2;
        Jxx = sum(Term1_Jxx + Term2_Jxx);
        
        Term1_Jyy = y_ue^2 ./ (r_n.^6);
        Term2_Jyy = (2*pi/lambda)^2 * (1./(r_n.^2) - 1./(r_1*r_n)).^2;
        Jyy = sum(Term1_Jyy + y_ue^2 * Term2_Jyy);
        
        Term1_Jxy = delta_n ./ (r_n.^6);
        Term2_Jxy = (2*pi/lambda)^2 * ((delta_n./(r_n.^2) - delta_1./(r_1*r_n)) .* (1./(r_n.^2) - 1./(r_1*r_n)));
        Jxy = sum( y_ue * (Term1_Jxy + Term2_Jxy) );
        
        % Calculo da potência total recebida para noise_power
        total_rx_power = 0;
        for n = 1:length(elements_x)
            d_n = sqrt( (x_ue-elements_x(n))^2 + y_ue^2 + elements_z(n)^2 );
            beta_n = (lambda/(4*pi*d_n))^2;
            total_rx_power = total_rx_power + P_tx*beta_n;
        end
        noise_power = total_rx_power / SNR;
        
        % Fator de normalização
        Pnorm = L * (lambda/(sqrt(8)*pi*sqrt(noise_power)))^2;
        Jxx_scaled = Pnorm * Jxx;
        Jyy_scaled = Pnorm * Jyy;
        Jxy_scaled = Pnorm * Jxy;
        
        % CRB em x e y
        CRB_den_x = Jxx_scaled - (Jxy_scaled^2 / Jyy_scaled);
        CRB_den_y = Jyy_scaled - (Jxy_scaled^2 / Jxx_scaled);
        CRB_x = 1 / CRB_den_x;
        CRB_y = 1 / CRB_den_y;
        CRB_norm = sqrt(CRB_x^2 + CRB_y^2);
        
        % Armazena os resultados no grid para a configuração atual
        CRB_norm_all(idx, i) = CRB_norm;
        CRB_x_all(idx, i)    = CRB_x;
        CRB_y_all(idx, i)    = CRB_y;
    end
end

%% Plot dos CRB versus altura para diferentes configurações (Nx = Ny)
colors = lines(N_configs);

% Combined CRB
figure; hold on;
for idx = 1:N_configs
    [min_val, min_idx] = min(CRB_norm_all(idx,:));
    z_min = z_vec(min_idx);
    plot(z_vec, CRB_norm_all(idx,:), 'LineWidth',2, 'Color', colors(idx,:), ...
        'DisplayName', ['Array = ' num2str(array_sizes(idx)) ' X ' num2str(array_sizes(idx))]);
end
xlabel('Antenna Mean Height (m)','Interpreter','tex');
ylabel('CRB');
grid on;
legend('Location','best');
set(gca, 'YScale', 'log');

elapsedTime = toc;
fprintf('Elapsed time: %.2f s\n', elapsedTime);
