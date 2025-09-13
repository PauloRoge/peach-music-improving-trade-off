%============================================================
%  PEACH-MUSIC – Comparação de variantes versus Elevação
%  (SNR fixo, varredura de elevação de 0 a 100 m)
%  Agora incluindo 4 variantes de PEACH com diferentes n_circ e n_hiper
%============================================================
clear; clc;
start_figure2;  % define Mx, Mz, d_x, d_z, lambda, L, MCS, deltaArea,...
                 % numIterNM, tol, x, y, n_hiper, n_circ, P_tx, alpha,...

% ---------- 1. Configurações --------------------------------------------
SNRdB_fixed    = 15;            
elev_vec       = 0:5:100;       
nelev          = numel(elev_vec);

% pré-alocação de RMSE
RMSE_peach     = zeros(1,nelev);
RMSE_peach1    = zeros(1,nelev);
RMSE_peach2    = zeros(1,nelev);
RMSE_peach3    = zeros(1,nelev);
RMSE_nm        = zeros(1,nelev);
RMSE_golden    = zeros(1,nelev);
CRBth          = zeros(1,nelev);

% pré-alocação de tempos (se quiser plotar também)
time_pch_avg   = zeros(1,nelev);
time_pch_avg1  = zeros(1,nelev);
time_pch_avg2  = zeros(1,nelev);
time_pch_avg3  = zeros(1,nelev);
time_nm_avg    = zeros(1,nelev);
time_gd_avg    = zeros(1,nelev);
time_un_avg    = zeros(1,nelev);

% ---------- 2. Loop principal (parfor em elevação) ----------------------
parfor k = 1:nelev
    elev_k = elev_vec(k);

    % --- recalcula sub-arrays para esta elevação -----------------------
    [URA, URA_x, URA_z, x_h, x_v, z_h, z_v] = ...
        subarrays(Mx, Mz, d_x, d_z, elev_k, lambda, 0);
    ref = URA(1,:);

    % acumuladores de erro e tempo
    err2_pc   = 0; err2_pc1 = 0; err2_pc2 = 0; err2_pc3 = 0;
    err2_nm   = 0; err2_gd  = 0; crb_sum  = 0;
    t_pc      = 0; t_pc1    = 0; t_pc2    = 0; t_pc3    = 0;
    t_nm      = 0; t_gd     = 0; t_un     = 0;

    for r = 1:MCS
        rng(r);                                
        pos = [ -50 + 100*rand , 13 + 37*rand , 0 ];

        % sinais recebidos
        [Yh, Yv, Y] = signals(pos, URA, lambda, L, ...
                              alpha, SNRdB_fixed, P_tx, Mx, Mz);

        % CRB
        crb_sum = crb_sum + crb(L, URA, pos, lambda, ...
                                P_tx, SNRdB_fixed, 2);

        % --- PEACH variante base (n_hiper, n_circ) ----------------------
        tic;
        [~, ~, est_pc] = peach_analitico( ...
            Yh, Yv, L, x, 24, x_h, z_h, x_v, z_v, ...
            ref, lambda, y, 12, pos);
        t_pc = t_pc + toc;

        % --- PEACH variante 1 (4× candidatos) ---------------------------
        tic;
        [~, ~, est_pc1] = peach_analitico( ...
            Yh, Yv, L, x, 98, x_h, z_h, x_v, z_v, ...
            ref, lambda, y, 46, pos);
        t_pc1 = t_pc1 + toc;

        % --- PEACH variante 2 (8× candidatos) ---------------------------
        tic;
        [~, ~, est_pc2] = peach_analitico( ...
            Yh, Yv, L, x, 192, x_h, z_h, x_v, z_v, ...
            ref, lambda, y, 98, pos);
        t_pc2 = t_pc2 + toc;

        % --- PEACH variante 3 (40× candidatos) --------------------------
        tic;
        [~, ~, est_pc3] = peach_analitico( ...
            Yh, Yv, L, x, 960, x_h, z_h, x_v, z_v, ...
            ref, lambda, y, 860, pos);
        t_pc3 = t_pc3 + toc;

        % --- Golden-PEACH -----------------------------------------------
        tic;
        [~, ~, est_gd] = peach_aurea( ...
            Yh, Yv, L, x, n_hiper, x_h, z_h, x_v, z_v, ...
            ref, lambda, y, n_circ, pos);
        t_gd = t_gd + toc;

        % --- subespaço de ruído (NM) -------------------------------------
        tic;
        Cov = (Y*Y')/L;
        [V,D] = eig(Cov); [~,idx] = sort(diag(D),'descend');
        Un = V(:, idx(2:end));  t_un = t_un + toc;

        % --- Nelder-Mead -----------------------------------------------
        tic;
        est_nm = nelder_mead(URA, est_pc1, Un, lambda, ref, ...
                             deltaArea, numIterNM, tol, false, x, y);
        t_nm = t_nm + toc;

        % acumula erros quadráticos
        err2_pc   = err2_pc   + norm(est_pc  - pos(1:2))^2;
        err2_pc1  = err2_pc1  + norm(est_pc1 - pos(1:2))^2;
        err2_pc2  = err2_pc2  + norm(est_pc2 - pos(1:2))^2;
        err2_pc3  = err2_pc3  + norm(est_pc3 - pos(1:2))^2;
        err2_gd   = err2_gd   + norm(est_gd  - pos(1:2))^2;
        err2_nm   = err2_nm   + norm(est_nm  - pos(1:2))^2;
    end

    % estatísticas por elevação
    RMSE_peach(k)  = sqrt(err2_pc  / MCS);
    RMSE_peach1(k) = sqrt(err2_pc1 / MCS);
    RMSE_peach2(k) = sqrt(err2_pc2 / MCS);
    RMSE_peach3(k) = sqrt(err2_pc3 / MCS);
    RMSE_golden(k) = sqrt(err2_gd  / MCS);
    RMSE_nm(k)     = sqrt(err2_nm  / MCS);
    CRBth(k)       = crb_sum    / MCS;

    time_pch_avg(k)  = t_pc  / MCS;
    time_pch_avg1(k) = t_pc1 / MCS;
    time_pch_avg2(k) = t_pc2 / MCS;
    time_pch_avg3(k) = t_pc3 / MCS;
    time_un_avg(k)   = t_un  / MCS;
    time_nm_avg(k)   = (t_nm  / MCS) + time_un_avg(k) + time_pch_avg(k);
    time_gd_avg(k)   = t_gd  / MCS;
end

% ---------- 3. Plot de RMSE vs Elevação -------------------------------
fig = figure('Units','centimeters','Position',[2 2 16 12]);
semilogy(elev_vec, RMSE_peach,  'x--', ...
         elev_vec, RMSE_peach1, 'x--', ...
         elev_vec, RMSE_peach2, 'x--', ...
         elev_vec, RMSE_peach3, 'x--', ...
         elev_vec, RMSE_nm,      'o-',  ...
         elev_vec, RMSE_golden,  's--', ...
         elev_vec, CRBth,        'k--', ...
         'LineWidth',1.5);
grid on;
xlabel('Elevação do arranjo (m)');
ylabel('RMSE (m)');
legend('PEACH 24–12','PEACH 96–48','PEACH 192–96','PEACH 960–480', ...
       'PEACH+NM','Golden-PEACH','CRB','Location','best');
