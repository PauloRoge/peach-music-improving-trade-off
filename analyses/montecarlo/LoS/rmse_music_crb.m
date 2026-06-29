% ==============================================================
% rmse_peach_cr.m
%  – compara MUSIC clássico com o CRB teórico (erro de norm corrigido)
% ==============================================================
clear; clc;
startup1;  % garante que as pastas com as funções estejam no path

% ------------------ 1. DEFINIÇÃO DE PARÂMETROS ------------------
% Dimensões da URA (linha × coluna)
Mx    = 8;           % número de elementos na dimensão x
Mz    = 8;           % número de elementos na dimensão z

% Espaçamento entre elementos (normalmente λ/2)
lambda = 0.1;        % comprimento de onda (exemplo: 0.1 m)
d_x    = lambda/2;   % espaçamento em x
d_z    = lambda/2;   % espaçamento em z

% Altura da URA (elevação constante z dos elementos)
elev   = 20;          % todos os elementos estão em z = elev

% Parâmetros do sinal e cenário
L      = 200;        % número de snapshots
P_tx   = 0.1;          % potência de transmissão (arbitrária)
alpha  = 1;          % parâmetro de atenuação (se houver)
MCS    = 10;       % número de execuções Monte‐Carlo

% Posição real do usuário [x; y; z]
pos    = [10; 30; 0]; % exemplo: usuário em x=5 m, y=10 m, z=0

% Definição da grade de busca (x e y)
x_min = -20;   x_max = 20;   n_x = 101;
y_min = -20;   y_max = 20;   n_y = 101;
x = linspace(x_min, x_max, n_x);   % vetor de candidatos em x
y = linspace(y_min, y_max, n_y);   % vetor de candidatos em y

% Parâmetros específicos para PEACH (não usados aqui, mas mantidos para contexto)
n_hiper = 100;   % número de pontos do grid na hipérbole (exemplo)
n_circ  = 100;   % número de pontos do círculo (exemplo)

% -----------------------------------------------------------------
% ------------------ 2. PREPARA ARRAY -----------------------------
% -----------------------------------------------------------------
[URA, URA_x, URA_z, x_h, x_v, z_h, z_v] = subarrays( ...
    Mx, Mz, d_x, d_z, elev, lambda, 0);
% URA é M×3, cada linha URA(k,:) = [x_k, y_k, z_k]

% Elemento de referência: primeiro elemento da URA
ref = URA(1,:);   % [x_ref, y_ref, z_ref]

% -----------------------------------------------------------------
% ------------------ 3. LAÇO DE SNRs E MONTE‐CARLO -----------------
% -----------------------------------------------------------------
SNR_dB_vec = 0:1:20;            % vetores de SNR em dB
RMSE      = zeros(size(SNR_dB_vec));
CRBth     = zeros(size(SNR_dB_vec));

for k = 1:numel(SNR_dB_vec)
    SNRdB = SNR_dB_vec(k);
    err2  = 0;   % soma dos erros quadráticos para este SNR
    
    % Cálculo do CRB teórico para este SNR
    CRBth(k) = crb(L, URA, pos, lambda, P_tx, SNRdB, alpha);
    
    % Monte‐Carlo
    for r = 1:MCS
        % Gera sinal (função signals)
        [Yh, Yv, Y] = signals(pos, URA, lambda, L, alpha, SNRdB, P_tx, Mx, Mz);
        
        % ------ ESTIMADOR MUSIC CLÁSSICO ------
        estMusic = music(Y, URA, L, lambda, ref, x, y);
        % estMusic é [x_est, y_est]
        
        % Soma do erro quadrático (apenas nas coordenadas x e y)
        % Correção: subtrair vetores linha 1×2 de vetores linha 1×2
        err2 = err2 + norm(estMusic - pos(1:2)')^2;
    end
    
    RMSE(k) = sqrt(err2 / MCS);
    fprintf('SNR = %+3d dB | RMSE = %.3f m | CRB = %.3f m\n', ...
            SNRdB, RMSE(k), CRBth(k));
end

% -----------------------------------------------------------------
% ------------------ 4. PLOT DOS RESULTADOS -----------------------
% -----------------------------------------------------------------
figure;
semilogy(SNR_dB_vec, RMSE, 'o-', 'LineWidth', 1.5); hold on;
semilogy(SNR_dB_vec, CRBth, 's--', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Erro (m)');
legend('RMSE Monte‐Carlo (MUSIC)', 'CRB teórico', 'Location', 'southwest');
title(sprintf('MUSIC clássico vs CRB  |  %dx%d URA,  L = %d,  MCS = %d', ...
              Mx, Mz, L, MCS));

% =================================================================
% ================ FUNÇÃO MUSIC CLÁSSICO ==========================
% =================================================================
function est = music(Y, URA, L, lambda, ref, x_grid, y_grid)
% MUSIC   Estima posição 2D [x, y] via MUSIC usando URA completa
%
%   est = music(Y, URA, L, lambda, ref, x_grid, y_grid) retorna um vetor [x_est, y_est]
%   que maximiza o pseudospectro de MUSIC.
%   - Y       : M×L matriz de snapshots completos (M = Mx×Mz)
%   - URA     : M×3, cada linha URA(k,:) = [x_k, y_k, z_k]
%   - L       : número de snapshots
%   - lambda  : comprimento de onda
%   - ref     : 1×3, [x_ref, y_ref, z_ref]
%   - x_grid  : vetor de candidatos em x (1×n_x)
%   - y_grid  : vetor de candidatos em y (1×n_y)
%
%   Retorna:
%   - est     : [x_est, y_est], estimativa de posição 2D

    % Número de elementos M
    [M, ~] = size(URA);
    
    % 1) Matriz de covariância
    R = (Y * Y') / L;
    
    % 2) Decomposição espectral de R
    [V, D] = eig(R);
    [~, ind] = sort(diag(D), 'descend');
    V = V(:, ind);  % ordena autovetores por autovalores decrescentes
    
    % 3) Subespaço‐ruído (assumindo 1 sinal)
    Un = V(:, 2:end);  % tamanho M×(M-1)
    
    % 4) Pré‐alocar matriz do pseudospectro
    n_x = length(x_grid);
    n_y = length(y_grid);
    Pmusic = zeros(n_x, n_y);
    
    x_ref = ref(1);
    y_ref = ref(2);
    z_ref = ref(3);
    
    % 5) Varredura sobre grid (x_i, y_j)
    for ix = 1:n_x
        x_i = x_grid(ix);
        for jy = 1:n_y
            y_j = y_grid(jy);
            
            % Altura do usuário (supondo plano z0 = 0)
            z0 = 0;
            
            % Distância do ponto candidato ao elemento de referência
            d_ref = sqrt((x_i - x_ref)^2 + (y_j - y_ref)^2 + (z0 - z_ref)^2);
            
            % Montar vetor de resposta a(x_i, y_j)
            a = zeros(M, 1);
            for k = 1:M
                xk = URA(k,1);
                yk = URA(k,2);
                zk = URA(k,3);
                d_k = sqrt((x_i - xk)^2 + (y_j - yk)^2 + (z0 - zk)^2);
                phi = - (2*pi / lambda) * (d_ref - d_k);
                a(k) = exp(1j * phi);
            end
            
            % 6) Avaliar pseudospectro de MUSIC
            denom = a' * (Un * Un') * a;
            Pmusic(ix, jy) = 1 / abs(denom);
        end
    end
    
    % 7) Busca pelo pico máximo
    [~, idx_max] = max(Pmusic(:));
    [ix_max, jy_max] = ind2sub(size(Pmusic), idx_max);
    
    % 8) Coordenada estimada
    est = [x_grid(ix_max), y_grid(jy_max)];
end
