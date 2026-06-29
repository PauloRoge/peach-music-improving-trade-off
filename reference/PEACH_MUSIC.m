clear ; close all; clc;
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Parâmetros e Geração do Cenário
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rng(0)
dd = 30;
% Usuário em (x, y, z=0)
% UEs = [-dd/sqrt(2), dd/sqrt(2), 0];

x_lim = 50;
y_lim = 50;
y_min = 10;

% % Sorteia posição do usuário
x_rand   = 30;
y_rand   = 30;
UEs = [x_rand, y_rand, 0];

% Sorteia posição do usuário de forma uniforme em um círculo de raio x_lim
% r = x_lim * sqrt(rand);     % raio (0 até x_lim)
% theta = 2*pi * rand;        % ângulo (0 até 2*pi)
% x_rand = r * cos(theta);
% y_rand = r * sin(theta);
% 
%UEs = [x_rand, y_rand, 0];


n_hiper = 816;
n_circ =  408;

% scale = 1.2;
% x_lim = scale * dd;
% y_lim = scale * dd;

% Dimensões do array 2D
Nx = 8; % horizontal
Ny = 8; % vertical
z_antenna = 20;
N = Nx * Ny;

L =10;           % número de time slots
SNR_dB = 5;
SNR = 10^(SNR_dB/10);

f_c = 15e9;
lambda = 3e8 / f_c;
d = lambda/2;

% Posição física das antenas 2D no plano xz (y=0)
x_offset = (Nx-1)*d/2;
z_offset = (Ny-1)*d/2;

elements_x_a = zeros(Nx, Ny);
elements_y_a = zeros(Nx, Ny);  % fixo em 0
elements_z_a = zeros(Nx, Ny);

for i = 1:Nx
    for j = 1:Ny
        elements_x_a(i,j) = (i-1)*d - x_offset;
        elements_y_a(i,j) = 0;  % y=0
        elements_z_a(i,j) = (j-1)*d - z_offset + z_antenna;
    end
end

figure;
scatter(elements_x_a(:), elements_z_a(:),'filled');
title('Antenas no plano XZ');
xlabel('x (m)'); ylabel('z (m)');
axis equal; grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) Geração do sinal: (2D), (horizontal) e (vertical)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X  = zeros(N, L);    % sinal completo (2D)
Xh = zeros(Nx, L);   % subarranjo horizontal  (linha j=1)
Xv = zeros(Ny, L);   % subarranjo vertical    (coluna i=1)

% Antena (1,1) como referência
x_ref = elements_x_a(1,1);
y_ref = elements_y_a(1,1);
z_ref = elements_z_a(1,1);

P_tx = 1;  % potência tx arbitrária

for t = 1:L
    sum_2d = zeros(N,1);
    sum_h  = zeros(Nx,1);
    sum_v  = zeros(Ny,1);
    total_rx_power = 0;

    UE = UEs(1,:);  % 1 usuário
    % Distância do usuário p/ antena ref:
    d_k1 = sqrt( (x_ref - UE(1))^2 + (y_ref - UE(2))^2 + (z_ref - UE(3))^2 );

    a_2d = zeros(N,1);
    a_h  = zeros(Nx,1);
    a_v  = zeros(Ny,1);

    for ix_ = 1:Nx
        for iy_ = 1:Ny
            idx = (ix_-1)*Ny + iy_;
            x_km = elements_x_a(ix_, iy_);
            y_km = elements_y_a(ix_, iy_);
            z_km = elements_z_a(ix_, iy_);

            d_km = sqrt((x_km - UE(1))^2 + (y_km - UE(2))^2 + (z_km - UE(3))^2);

            beta_ = (lambda/(4*pi*d_km))^2;
            total_rx_power = total_rx_power + beta_;

            phase_ = -(2*pi/lambda)*(d_k1 - d_km);
            a_2d(idx) = sqrt(beta_)*exp(1j*phase_);

            % Linha j=1 => subarranjo horizontal
            if (iy_==1)
                a_h(ix_) = sqrt(beta_)*exp(1j*phase_);
            end
            % Coluna i=1 => subarranjo vertical
            if (ix_==1)
                a_v(iy_) = sqrt(beta_)*exp(1j*phase_);
            end
        end
    end

    s_ = (randn + 1j*randn);
    sum_2d = sum_2d + a_2d*s_;
    sum_h  = sum_h  + a_h*s_;
    sum_v  = sum_v  + a_v*s_;

    noise_power = total_rx_power / SNR;

    X(:,t)  = sum_2d + sqrt(noise_power/2)*(randn(N,1)+1j*randn(N,1));
    Xh(:,t) = sum_h  + sqrt(noise_power/2)*(randn(Nx,1)+1j*randn(Nx,1));
    Xv(:,t) = sum_v  + sqrt(noise_power/2)*(randn(Ny,1)+1j*randn(Ny,1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3) Matrizes de covariância e subespaço de ruído
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rxx  = (1/L)*(X  * X' );
Rxxh = (1/L)*(Xh * Xh');
Rxxv = (1/L)*(Xv * Xv');

[U,  ~,~] = svd(Rxx );
[Uh, ~,~] = svd(Rxxh);
[Uv, ~,~] = svd(Rxxv);

% Sabemos K=1 => subespaço de ruído
Un  = U(:,  2:end);
Unh = Uh(:, 2:end);
Unv = Uv(:, 2:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PEACH-MUSIC:
%% 4) Subarranjo horizontal => HIPÉRBOLE => extrair x_peak e depois Delta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_candidates = linspace(-x_lim, x_lim, n_hiper);
Pmusic_h = zeros(size(x_candidates));

for i = 1:length(x_candidates)
    x_cand = x_candidates(i);
    a_h_ = zeros(Nx,1);

    % Distância p/ antena ref (y=0)
    d_ref_ = sqrt( (x_ref - x_cand)^2 + z_ref^2 );

    for ix_ = 1:Nx
        x_km = elements_x_a(ix_,1);
        z_km = elements_z_a(ix_,1);
        d_km_ = sqrt( (x_km - x_cand)^2 + z_km^2 );

        phase_ = -(2*pi/lambda)*( d_ref_ - d_km_ );
        a_h_(ix_) = exp(1j*phase_);
    end
    Pmusic_h(i) = 1 / abs( a_h_'*(Unh*Unh')*a_h_ );
end

[~, idx_max_h] = max(Pmusic_h);
x_peak_h = x_candidates(idx_max_h);

% Calculamos \Delta = dF1 - dF2 com base nas antenas extremas
F1x = elements_x_a(1,1);
F1z = elements_z_a(1,1);
F2x = elements_x_a(Nx,1);
F2z = elements_z_a(Nx,1);

dF1 = sqrt( (x_peak_h - F1x)^2 + (0 - 0)^2 + (0 - F1z)^2 );
dF2 = sqrt( (x_peak_h - F2x)^2 + (0 - 0)^2 + (0 - F2z)^2 );
Delta_est = dF1 - dF2;

fprintf('\n[Hiperbole] x_peak=%.2f => Delta=%.2f\n', x_peak_h, Delta_est);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5) Subarranjo vertical => CÍRCULO => extrair y_peak => R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y_candidates = linspace(0, y_lim, n_circ);
Pmusic_v = zeros(size(y_candidates));

for i = 1:length(y_candidates)
    y_cand = y_candidates(i);
    a_v_ = zeros(Ny,1);

    d_ref_ = sqrt( (x_ref - 0)^2 + (y_ref - y_cand)^2 + (z_ref - 0)^2 );

    for iy_ = 1:Ny
        x_km = elements_x_a(1,iy_);
        z_km = elements_z_a(1,iy_);
        d_km_ = sqrt( (x_km - 0)^2 + (0 - y_cand)^2 + (z_km - 0)^2 );

        phase_ = -(2*pi/lambda)*( d_ref_ - d_km_ );
        a_v_(iy_) = exp(1j*phase_);
    end
    Pmusic_v(i) = 1 / abs( a_v_'*(Unv*Unv')*a_v_ );
end

[~, idx_max_v] = max(Pmusic_v);
y_peak_v = y_candidates(idx_max_v);

R_est = y_peak_v;  % "raio"

fprintf('[Circulo] y_peak=%.2f => R=%.2f\n', y_peak_v, R_est);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6) Resolver interseção: HIPÉRBOLE vs CÍRCULO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% =============== MÉTODO NUMÉRICO (solve) ===============
syms xs ys real

expr_d1 = sqrt( (xs - F1x)^2 + (ys - 0)^2 + (0 - F1z)^2 );
expr_d2 = sqrt( (xs - F2x)^2 + (ys - 0)^2 + (0 - F2z)^2 );
eq_hip  = expr_d1 - expr_d2 - Delta_est;      % Diferença de distâncias 3D
eq_circ = xs^2 + ys^2 - R_est^2;             % Círculo no plano (x,y)

sol = solve([eq_hip, eq_circ],[xs, ys],'Real',true);
x_sol = double(sol.xs);
y_sol = double(sol.ys);

fprintf('\nSoluções da Interseção (Hipérbole & Círculo) [solve]:\n');
for i=1:length(x_sol)
    fprintf('  (%.2f, %.2f)\n', x_sol(i), y_sol(i));
end

% =============== MÉTODO ANALÍTICO ===============
% Supondo que:
%   - F1 = (-c, 0, zA) e F2 = (+c, 0, zA)
%   - Delta_est = d1 - d2
%   - R_est = d_xy
%   - c_analitico = (F2x - F1x)/2
%   - zA = F1z (assumindo F1z = F2z)

c_analitico = (F2x - F1x)/2;   % deve ser >0 se F2x > F1x
zA          = F1z;           % antenas na mesma altura

% Fórmula: x^2 = (Delta^2 / (16*c^2)) * [4(c^2 + zA^2 + R_est^2) - Delta^2]
temp_   = 4*(c_analitico^2 + zA^2 + R_est^2) - Delta_est^2;
x2_anal = (Delta_est^2/(16 * c_analitico^2)) * temp_;

x_sol2 = [];  % Vetores de solução analítica
y_sol2 = [];

if x2_anal < 0
    fprintf('\n[Solução Analítica] x^2 < 0 => sem interseção real.\n');
else
    % x => ± sqrt(x2_anal)
    x_pos = +sqrt(x2_anal);
    x_neg = -sqrt(x2_anal);

    % y^2 = R_est^2 - x^2
    y2_from_xpos = R_est^2 - x_pos^2;
    y2_from_xneg = R_est^2 - x_neg^2;

    % Combinações possíveis de sinais
    cand_x = [x_pos, x_pos, x_neg, x_neg];
    cand_y = [ ...
        +sqrt(max(y2_from_xpos,0)), -sqrt(max(y2_from_xpos,0)), ...
        +sqrt(max(y2_from_xneg,0)), -sqrt(max(y2_from_xneg,0)) ];

    % Filtrar NaN/Inf
    valid_idx = find(~isnan(cand_y) & ~isinf(cand_y));
    x_sol2 = cand_x(valid_idx);
    y_sol2 = cand_y(valid_idx);
end

fprintf('\nSoluções da Interseção (Hipérbole & Círculo) [Analítico]:\n');
for ii = 1:length(x_sol2)
    fprintf('  (%.2f, %.2f)\n', x_sol2(ii), y_sol2(ii));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7) Escolher a solução mais próxima da real
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(x_sol)
    fprintf('\nNão foram encontradas soluções reais para a interseção.\n');
    pos_est = [NaN, NaN];
    best_dist = Inf;
else
    % Se há pelo menos 1 solução, faça a busca do "best_idx"
    UE_real = [UEs(1), UEs(2)];
    best_idx = 1;
    best_dist = Inf;
    for i = 1:length(x_sol)
        cand = [x_sol(i), y_sol(i)];
        dist_ = norm(cand - UE_real);
        if dist_ < best_dist
            best_dist = dist_;
            best_idx = i;
        end
    end
    
    pos_est = [x_sol(best_idx), y_sol(best_idx)];
    fprintf('\nPosição REAL do usuário: (%.2f, %.2f)\n', UE_real(1), UE_real(2));
    fprintf('Posição ESTIMADA: (%.2f, %.2f)\n', pos_est(1), pos_est(2));
    fprintf('Erro euclidiano = %.2f (m)\n', best_dist);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (Trecho Adicional) Geração dos Mapas 2D de Pseudo-Espectro
%%   1) Subarranjo Horizontal
%%   2) Subarranjo Vertical
%%   3) Soma dos dois
%%
%% A suposição é que:
%%  - Nx, Ny, lambda, elements_x_a, elements_z_a, etc. já existem
%%  - Unh, Unv já calculados
%%  - Você deseja varrer x e y em uma grade 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Defina a resolução e a faixa de varredura:
Nx_plot = 100;  
Ny_plot = 100;

x_vec = linspace(-x_lim, x_lim, Nx_plot);
y_vec = linspace(0, y_lim, Ny_plot);  % <-- agora só y>=0

Pmap_h = zeros(Ny_plot, Nx_plot);  
Pmap_v = zeros(Ny_plot, Nx_plot);

% -- Subarranjo Horizontal --
for ix = 1:Nx_plot
    for iy = 1:Ny_plot
        x_cand = x_vec(ix);
        y_cand = y_vec(iy);

        a_h_ = zeros(Nx,1);
        d_ref_ = sqrt( (x_ref - x_cand)^2 + (y_ref - y_cand)^2 + z_ref^2 );
        
        for n_ = 1:Nx
            x_ant = elements_x_a(n_,1);
            z_ant = elements_z_a(n_,1);
            % y=0 para antenas horizontais
            d_km_ = sqrt( (x_ant - x_cand)^2 + (0 - y_cand)^2 + z_ant^2 );
            phase__ = -(2*pi/lambda)*( d_ref_ - d_km_ );
            a_h_(n_) = exp(1j*phase__);
        end
        Pmap_h(iy, ix) = 1 / abs( a_h_'*(Unh*Unh')*a_h_ );
    end
end

% -- Subarranjo Vertical --
for ix = 1:Nx_plot
    for iy = 1:Ny_plot
        x_cand = x_vec(ix);
        y_cand = y_vec(iy);

        a_v_ = zeros(Ny,1);
        d_ref_ = sqrt( (x_ref - x_cand)^2 + (y_ref - y_cand)^2 + (z_ref - 0)^2 );

        for m_ = 1:Ny
            x_ant = elements_x_a(1,m_);
            z_ant = elements_z_a(1,m_);
            % x=0 para antenas verticais
            d_km_ = sqrt( (x_ant - x_cand)^2 + (0 - y_cand)^2 + (z_ant - 0)^2 );
            phase__ = -(2*pi/lambda)*( d_ref_ - d_km_ );
            a_v_(m_) = exp(1j*phase__);
        end
        Pmap_v(iy, ix) = 1 / abs( a_v_'*(Unv*Unv')*a_v_ );
    end
end

Pmap_sum = Pmap_h + Pmap_v;
Pmap_mul = Pmap_h .* Pmap_v;

figure('Name','Mapas de Calor - Somente y >= 0','Color',[1 1 1]);

% Subplot 1
subplot(2,2,1);
imagesc(x_vec, y_vec, 10*log10(Pmap_h)/ max(max( 10*log10(Pmap_h))));
set(gca,'YDir','normal'); colorbar;
title('Horizontal (Hipérbole)');
xlabel('x (m)'); ylabel('y (m)');
hold on; plot(UEs(1), UEs(2), 'ro','MarkerFaceColor','r','MarkerSize',7);
hold off;

% Subplot 2
subplot(2,2,2);
imagesc(x_vec, y_vec, 10*log10(Pmap_v)/max(max( 10*log10(Pmap_v))));
set(gca,'YDir','normal'); colorbar;
title('Vertical (Círculo)');
xlabel('x (m)'); ylabel('y (m)');
hold on; plot(UEs(1), UEs(2), 'ro','MarkerFaceColor','r','MarkerSize',7);
hold off;

% Subplot 3
subplot(2,2,3);
imagesc(x_vec, y_vec, 10*log10(Pmap_sum)/max(max( 10*log10(Pmap_sum))));
set(gca,'YDir','normal'); colorbar;
title('Soma (Hipérbole + Círculo)');
xlabel('x (m)'); ylabel('y (m)');
hold on; plot(UEs(1), UEs(2), 'ro','MarkerFaceColor','r','MarkerSize',7);
hold off;

% Subplot 3
subplot(2,2,4);
imagesc(x_vec, y_vec, 10*log10(Pmap_mul)/max(max( 10*log10(Pmap_mul))));
set(gca,'YDir','normal'); colorbar;
title('Soma (Hipérbole + Círculo)');
xlabel('x (m)'); ylabel('y (m)');
hold on; plot(UEs(1), UEs(2), 'ro','MarkerFaceColor','r','MarkerSize',7);
hold off;



% Subplot 1
figure(1337)
imagesc(x_vec, y_vec, 10*log10(Pmap_h)/ max(max( 10*log10(Pmap_h))));
set(gca,'YDir','normal'); colorbar;
xlabel('x (m)'); ylabel('y (m)');
hold on; plot(UEs(1), UEs(2), 'ro','MarkerFaceColor','r','MarkerSize',7);
hold off;

% Subplot 2
figure(1338)
imagesc(x_vec, y_vec, 10*log10(Pmap_v)/max(max( 10*log10(Pmap_v))));
set(gca,'YDir','normal'); colorbar;
xlabel('x (m)'); ylabel('y (m)');
hold on; plot(UEs(1), UEs(2), 'ro','MarkerFaceColor','r','MarkerSize',7);
hold off;

% Subplot 3
figure(1339)
imagesc(x_vec, y_vec, 10*log10(Pmap_sum)/max(max( 10*log10(Pmap_sum))));
set(gca,'YDir','normal'); colorbar;
xlabel('x (m)'); ylabel('y (m)');
hold on; plot(UEs(1), UEs(2), 'ro','MarkerFaceColor','r','MarkerSize',7);
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pseudo-espectro 2D (todas as antenas) e mapa de calor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vamos usar a mesma resolução Nx_plot, Ny_plot e os mesmos x_vec, y_vec
% já definidos anteriormente (ex.: Nx_plot=100, Ny_plot=100, etc.).
% Se quiser outra faixa em y (inclusive valores negativos), basta ajustar y_vec.

Pmap_2d = zeros(Ny_plot, Nx_plot);

for ix_ = 1:Nx_plot
    for iy_ = 1:Ny_plot
        x_cand = x_vec(ix_);
        y_cand = y_vec(iy_);

        % Monta o steering vector para TODAS as antenas (Nx x Ny)
        a_2d_ = zeros(N,1);
        d_ref_ = sqrt( (x_ref - x_cand)^2 + (y_ref - y_cand)^2 + z_ref^2 );

        idx = 1;
        for ixx = 1:Nx
            for iyy = 1:Ny
                dx_ = elements_x_a(ixx, iyy);
                dy_ = elements_y_a(ixx, iyy);
                dz_ = elements_z_a(ixx, iyy);

                % Distância do candidato (x_cand, y_cand) até a antena (dx_, dy_, dz_)
                d_km_ = sqrt( (dx_ - x_cand)^2 + (dy_ - y_cand)^2 + (dz_ - 0)^2 );

                % Fase relativa à antena de referência
                ph_ = -(2*pi/lambda)*( d_ref_ - d_km_ );

                a_2d_(idx) = exp(1j*ph_);
                idx = idx + 1;
            end
        end

        % Pseudo-espectro MUSIC: 1 / |a^H * Un * Un^H * a|
        % Obs: "Un" vem da SVD de Rxx = (1/L)*X*X', referente ao array completo
        Pmap_2d(iy_, ix_) = 1 / abs( a_2d_'*(Un*Un')*a_2d_ );
    end
end

% Faz o plot em dB (normalizado para o máximo)
figure('Name','Pseudo-espectro 2D (todas as antenas)','Color',[1 1 1]);
imagesc(x_vec, y_vec, 10*log10(Pmap_2d) - 10*log10(max(Pmap_2d(:))));
set(gca,'YDir','normal'); 
colorbar;
title('Pseudo-espectro - URA Completa (Nx x Ny)');
xlabel('x (m)'); 
ylabel('y (m)');

hold on; 
    % 1) Posição REAL do usuário (vermelho)
    plot(UEs(1), UEs(2), 'ro','MarkerFaceColor','r','MarkerSize',7); 
    
    % 2) Encontra índice do valor máximo em Pmap_2d (pico do pseudo-espectro)
    [~, idx_max] = max(Pmap_2d(:)); 
    [iy_max, ix_max] = ind2sub(size(Pmap_2d), idx_max); 
    x_max = x_vec(ix_max);
    y_max = y_vec(iy_max);
    
    % 3) Marca o pico do pseudo-espectro (estrela preta)
    plot(x_max, y_max, 'k*', 'MarkerSize',10, 'LineWidth',1.5);
    
    % 4) Posição ESTIMADA do PEACH (verde) - use 'pos_est' já calculado
    plot(pos_est(1), pos_est(2), 'go', 'MarkerFaceColor','g','MarkerSize',7);
hold off;






toc;
