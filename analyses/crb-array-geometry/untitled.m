clear; close all; clc;

%% Configurações de arranjos desejadas
configs = [ ...
    8 8;
   16 4;
   4 16;
   2 32;
   32 2;
];

%% Região de busca
x_min = -100; x_max = 100;
y_min = 10;  y_max = 110;
x_vec = linspace(x_min, x_max, 200);
y_vec = linspace(y_min, y_max, 200);
[X_grid, Y_grid] = meshgrid(x_vec, y_vec);

%% Parâmetros fixos
f_c = 15e9;
lambda = 3e8 / f_c;
d = 0.5 * lambda;
z_antenna = 20;
P_tx = 1;
SNR_dB = 15;
SNR = 10^(SNR_dB / 10);
L = 10;

%% Loop sobre configurações
for cfg = 1:size(configs,1)
    Nx = configs(cfg,1);
    Ny = configs(cfg,2);
    N = Nx * Ny;
    
    % Posicionamento dos elementos
    elements_x = zeros(N,1);
    elements_z = zeros(N,1);
    x_offset = ((Nx - 1) * d / 2);
    z_offset = ((Ny - 1) * d / 2);
    idx = 1;
    for i = 1:Nx
        for j = 1:Ny
            elements_x(idx) = (i - 1) * d - x_offset;
            elements_z(idx) = (j - 1) * d - z_offset + z_antenna;
            idx = idx + 1;
        end
    end

    % Normalização fixa como no CRB_HEAT
    noise_power = P_tx / SNR;
    Pnorm = L * (lambda / (sqrt(8) * pi * sqrt(noise_power)))^2;

    % Inicializa matrizes de CRB
    CRB_x_map = zeros(size(X_grid));
    CRB_y_map = zeros(size(X_grid));
    CRB_tot = zeros(size(X_grid));

    % Loop sobre a malha (x, y)
    for i = 1:length(y_vec)
        y_ue = y_vec(i);
        for j = 1:length(x_vec)
            x_ue = x_vec(j);

            x_1 = elements_x(1);
            z_1 = elements_z(1);
            r_1 = sqrt((x_ue - x_1)^2 + y_ue^2 + z_1^2);
            delta_1 = (x_ue - x_1);

            delta_n = (x_ue - elements_x);
            r_n = sqrt((x_ue - elements_x).^2 + y_ue^2 + elements_z.^2);

            Term1_Jxx = (delta_n.^2) ./ r_n.^6;
            Term2_Jxx = (2 * pi / lambda)^2 * (delta_n ./ r_n.^2 - delta_1 ./ (r_1 .* r_n)).^2;
            Jxx = sum(Term1_Jxx + Term2_Jxx);

            Term1_Jyy = y_ue^2 ./ r_n.^6;
            Term2_Jyy = (2 * pi / lambda)^2 * (1 ./ r_n.^2 - 1 ./ (r_1 .* r_n)).^2;
            Jyy = sum(Term1_Jyy + y_ue^2 * Term2_Jyy);

            Term1_Jxy = delta_n ./ r_n.^6;
            Term2_Jxy = (2 * pi / lambda)^2 * ...
                (delta_n ./ r_n.^2 - delta_1 ./ (r_1 .* r_n)) .* ...
                (1 ./ r_n.^2 - 1 ./ (r_1 .* r_n));
            Jxy = sum(y_ue * (Term1_Jxy + Term2_Jxy));

            Jxx_scaled = Pnorm * Jxx;
            Jyy_scaled = Pnorm * Jyy;
            Jxy_scaled = Pnorm * Jxy;

            CRB_den_x = Jxx_scaled - (Jxy_scaled^2 / Jyy_scaled);
            CRB_den_y = Jyy_scaled - (Jxy_scaled^2 / Jxx_scaled);

            if CRB_den_x > 0
                CRB_x_map(i,j) = 1 / CRB_den_x;
            else
                CRB_x_map(i,j) = NaN;
            end

            if CRB_den_y > 0
                CRB_y_map(i,j) = 1 / CRB_den_y;
            else
                CRB_y_map(i,j) = NaN;
            end

            if CRB_den_x > 0 && CRB_den_y > 0
                CRB_tot(i,j) = sqrt(CRB_x_map(i,j)^2 + CRB_y_map(i,j)^2);
            else
                CRB_tot(i,j) = NaN;
            end
        end
    end

    %% Plots
    % Mapa binário CRB_x > CRB_y
    CRB_diff = CRB_x_map - CRB_y_map;
    bin_map = zeros(size(CRB_diff));
    bin_map(CRB_diff > 0) = 1;

    figure;
    imagesc(x_vec, y_vec, bin_map);
    set(gca,'YDir','normal');
    xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('CRB x > y - %dx%d', Nx, Ny));
    colorbar;

    % CRB total
    figure;
    imagesc(x_vec, y_vec, log10(CRB_tot));
    clim([3 10]); %Evita interpretações enganosas (“parece melhor, mas só mudou a escala”).
    set(gca,'YDir','normal');
    xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('CRB total - %dx%d', Nx, Ny));
    colorbar;

    % CRB y
    figure;
    imagesc(x_vec, y_vec, log10(CRB_y_map));
    clim([3 10]);  %cor azul é sempre log10(CRB) ≈ 3 e amarela é sempre ≈ 10
    set(gca,'YDir','normal');
    xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('CRB y - %dx%d', Nx, Ny));
    colorbar;

    % CRB x
    figure;
    imagesc(x_vec, y_vec, log10(CRB_x_map));
    clim([3 10]);  % igual ao CRB_HEAT.m
    set(gca,'YDir','normal');
    xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('CRB x - %dx%d', Nx, Ny));
    colorbar;
end
