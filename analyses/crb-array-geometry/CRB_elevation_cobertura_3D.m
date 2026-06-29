clear; close all; clc;
tic;

%% Parameters
dd      = 50;
scale   = 1.5;
L       = 10;
x_min   = -100; x_max = 100;      % User x range
y_min   = 120;  y_max = 200;       % User y range


%nearfield
% x_min   = -100; x_max = 100;      % User x range
% y_min   = 120;  y_max = 200;       % User y range

%transição
% x_min   = -100; x_max = 100;      % User x range
% y_min   = 120;  y_max = 200;       % User y range

% farfield
% x_min   = -50; x_max = 50;      % User x range
% y_min   = 120;  y_max = 200;       % User y range
Nsub    = 1;
dspace  = 0.5;
f_c     = 15e9;
lambda  = 3e8 / f_c;
d       = 0.5 * lambda;
P_tx    = 1;
SNR_dB  = 30;
SNR     = 10^(SNR_dB/10);

% Grid for user positions
nPts = 32;
x_vec = linspace(x_min, x_max, nPts);
y_vec = linspace(y_min, y_max, nPts);
[Xgrid, Ygrid] = meshgrid(x_vec, y_vec);

% Heights to simulate
z_values = [0 10 20 30 40 50];
num_z = length(z_values);

% Configurações de array: [Nx Ny]
array_configs = [32 2; 16 4; 8 8; 4 16; 2 32];

%% Loop sobre configurações de array
for cfg = 1:size(array_configs,1)
    Nx = array_configs(cfg,1);
    Ny = array_configs(cfg,2);
    
    figure; hold on;
    colors = parula(num_z);
    title(sprintf('CRB for Array %dx%d', Nx, Ny));

    %% Loop sobre alturas
    for idx = 1:num_z
        z_ant = z_values(idx);
        CRB_surf = zeros(size(Xgrid));
        
        for i = 1:nPts
            for j = 1:nPts
                x_ue = Xgrid(i,j);
                y_ue = Ygrid(i,j);
                
                % Gerar posições dos elementos da array
                elements_x_a = zeros(Nsub*Nx, Ny);
                elements_y_a = zeros(Nsub*Nx, Ny);
                elements_z_a = zeros(Nsub*Nx, Ny);
                x_offset = ((Nx-1)*d/2) + (Nsub-1)*(Nx*d+dspace)/2;
                z_offset = ((Ny-1)*d/2);
                for k_sub = 1:Nsub
                    for i_ant = 1:Nx
                        for j_ant = 1:Ny
                            idx_elem = (k_sub-1)*Nx + i_ant;
                            elements_x_a(idx_elem,j_ant) = (i_ant-1)*d + (k_sub-1)*(Nx*d+dspace) - x_offset;
                            elements_y_a(idx_elem,j_ant) = 0;
                            elements_z_a(idx_elem,j_ant) = (j_ant-1)*d - z_offset + z_ant;
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

                % FIM
                delta_1 = abs(x_ue - x_k1);
                r_1 = sqrt( (x_ue - x_k1)^2 + y_ue^2 + z_k1^2 );
                delta_n = abs(x_ue - elements_x);
                r_n = sqrt( (x_ue - elements_x).^2 + y_ue^2 + elements_z.^2 );

                Term1_Jxx = (delta_n.^2) ./ (r_n.^6);
                Term2_Jxx = (2*pi/lambda)^2 * ( delta_n./(r_n.^2) - delta_1./(r_1*r_n) ).^2;
                Jxx = sum(Term1_Jxx + Term2_Jxx);

                Term1_Jyy = y_ue^2 ./ (r_n.^6);
                Term2_Jyy = (2*pi/lambda)^2 * ( 1./(r_n.^2) - 1./(r_1*r_n) ).^2;
                Jyy = sum(Term1_Jyy + y_ue^2 * Term2_Jyy);

                Term1_Jxy = delta_n ./ (r_n.^6);
                Term2_Jxy = (2*pi/lambda)^2 * ((delta_n./(r_n.^2) - delta_1./(r_1*r_n)) .* (1./(r_n.^2) - 1./(r_1*r_n)));
                Jxy = sum( y_ue * (Term1_Jxy + Term2_Jxy) );

                % Potência recebida
                total_rx_power = 0;
                for n = 1:length(elements_x)
                    d_n = sqrt((x_ue-elements_x(n))^2 + y_ue^2 + elements_z(n)^2);
                    beta_n = (lambda/(4*pi*d_n))^2;
                    total_rx_power = total_rx_power + P_tx * beta_n;
                end
                noise_power = total_rx_power / SNR;

                % Normalização
                Pnorm = L * (lambda/(sqrt(8)*pi*sqrt(noise_power)))^2;
                Jxx_scaled = Pnorm * Jxx;
                Jyy_scaled = Pnorm * Jyy;
                Jxy_scaled = Pnorm * Jxy;

                % CRB combinado
                CRB_den_x = Jxx_scaled - (Jxy_scaled^2 / Jyy_scaled);
                CRB_den_y = Jyy_scaled - (Jxy_scaled^2 / Jxx_scaled);
                CRB_x = 1/CRB_den_x;
                CRB_y = 1/CRB_den_y;
                CRB_norm = sqrt(CRB_x^2 + CRB_y^2);
                
                CRB_surf(i,j) = CRB_norm;
            end
        end

        % Plot da superfície
        h = surf(Xgrid, Ygrid, CRB_surf, 'EdgeColor','none');
        h.FaceAlpha = 0.7;
        h.FaceColor = colors(idx,:);
        h.DisplayName = ['z = ' num2str(z_ant) ' m'];
    end

    xlabel('User x (m)');
    ylabel('User y (m)');
    zlabel('Combined CRB');
    legend('Location','best');
    colorbar;
    grid on;
    view(45,30);
    set(gca, 'ZScale', 'log');
    %zticks(10.^(-2:1:4));
    %zlim([1e-1 1e1]);
end

elapsedTime = toc;
fprintf('Elapsed time: %.2f s\n', elapsedTime);
