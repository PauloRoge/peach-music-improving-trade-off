%============================================================
%  PEACH-MUSIC – Comparação de variantes versus Elevação
%  (SNR fixo, varredura de elevação de 0 a 100 m)
%  Baseado em figure2.m
%============================================================
clear; clc;
start_figure2;                               % parâmetros gerais do projeto

% ---------- 1. Configurações --------------------------------------------
SNRdB_fixed = 15;            % SNR constante (dB)
elev_vec    = 0:5:100;       % elevação do array (m) — passo de 5 m
nelev       = numel(elev_vec);

RMSE_peach  = zeros(1,nelev);
RMSE_nm     = zeros(1,nelev);
% RMSE_hj     = zeros(1,nelev);
% RMSE_sbplx  = zeros(1,nelev);
RMSE_golden = zeros(1,nelev);
CRBth       = zeros(1,nelev);

time_pch_avg = zeros(1,nelev);
time_nm_avg  = zeros(1,nelev);
% time_hj_avg  = zeros(1,nelev);
% time_sb_avg  = zeros(1,nelev);
time_gd_avg  = zeros(1,nelev);

time_un_avg  = zeros(1,nelev);

% ---------- 2. Loop principal (parfor em elevação) ----------------------
parfor k = 1:nelev
    elev_k = elev_vec(k);

    % --- sub-arrays para esta elevação ----------------------------------
    [URA, URA_x, URA_z, x_h, x_v, z_h, z_v] = ...
        subarrays(Mx, Mz, d_x, d_z, elev_k, lambda, 0);
    ref = URA(1,:);

    % acumuladores de erro e tempo
    err2_pc = 0; err2_nm = 0; err2_hj = 0; err2_sb = 0; err2_gd = 0;
    t_pc = 0; t_nm = 0; 
    %t_hj = 0; t_sb = 0;
    t_gd = 0; t_un = 0;
    crb_sum = 0;

    for r = 1:MCS
        rng(r);                                    % reprodutibilidade
        pos = [ -50 + 100*rand , 13 + 37*rand , 0 ];

        % ---------- Sinais ----------------------------------------------
        [Yh, Yv, Y] = signals_los(pos, URA, lambda, L, ...
                              alpha, SNRdB_fixed, P_tx, Mx, Mz);

        % ---------- CRB --------------------------------------------------
        crb_sum = crb_sum + crb(L, URA, pos, lambda, ...
                                P_tx, SNRdB_fixed, 2);

        % ---------- PEACH (coarse) --------------------------------------
        tic;
        [~, ~, est_peach] = peach( ...
            Yh, Yv, L, x, n_hiper, x_h, z_h, x_v, z_v, ...
            ref, lambda, y, n_circ, pos);
        t_pc = t_pc + toc;

        % ---------- Golden-PEACH ----------------------------------------
        tic;
        [~, ~, est_golden] = peach_golden( ...
            Yh, Yv, L, x, n_hiper, x_h, z_h, x_v, z_v, ...
            ref, lambda, y, n_circ, pos);
        t_gd = t_gd + toc;

        % ---------- Sub-espaço de ruído ---------------------------------
        tic;
        Cov = (Y*Y')/L;
        [V,D] = eig(Cov); [~,idx] = sort(diag(D),'descend');
        Un = V(:, idx(2:end));                       % fonte única
        t_un = t_un + toc;

        % ---------- Hooke-Jeeves ----------------------------------------
        % tic;
        % hj_est = hooke_jeeves(URA, est_peach, Un, lambda, ...
        %                       ref, deltaArea, 10, tol, x, y, true);
        % t_hj = t_hj + toc;

        % ---------- Nelder-Mead -----------------------------------------
        tic;
        est_nm = nelder_mead(URA, est_peach, Un, lambda, ref, ...
                             deltaArea, numIterNM, tol, false, x, y);
        t_nm = t_nm + toc;

        % ---------- Subplex ---------------------------------------------
        % try
        %     tic;
        %     est_sb = subplex_wrapper(URA, est_peach, Un, lambda, ref, ...
        %                              x, y, tol_subplex, 2);
        %     t_sb = t_sb + toc;
        % catch
        %     est_sb = est_nm;                         % fallback
        % end

        % ---------- Acumula erros ---------------------------------------
        err2_pc = err2_pc + norm(est_peach  - pos(1:2))^2;
        err2_nm = err2_nm + norm(est_nm     - pos(1:2))^2;
        % err2_hj = err2_hj + norm(hj_est     - pos(1:2))^2;
        % err2_sb = err2_sb + norm(est_sb     - pos(1:2))^2;
        err2_gd = err2_gd + norm(est_golden - pos(1:2))^2;
    end

    % ---------- Estatísticas por elevação -------------------------------
    RMSE_peach(k)  = sqrt(err2_pc / MCS);
    RMSE_nm(k)     = sqrt(err2_nm / MCS);
    % RMSE_hj(k)     = sqrt(err2_hj / MCS);
    % RMSE_sbplx(k)  = sqrt(err2_sb / MCS);
    RMSE_golden(k) = sqrt(err2_gd / MCS);
    CRBth(k)       = crb_sum / MCS;

    time_un_avg(k)  = t_un / MCS;
    time_pch_avg(k) = t_pc / MCS;
    time_nm_avg(k)  = (t_nm / MCS) + time_un_avg(k) + time_pch_avg(k);
    % time_hj_avg(k)  = (t_hj / MCS) + time_un_avg(k) + time_pch_avg(k);
    % time_sb_avg(k)  = (t_sb / MCS) + time_un_avg(k) + time_pch_avg(k);
    time_gd_avg(k)  = t_gd / MCS;

    fprintf(['Elev %3.0f m | RMSE→ PCH %.2f  NM %.2f  HJ %.2f  SB %.2f  ' ...
             'GD %.2f | CRB %.2f  (avg de %.0f exec.)\n'], ...
            elev_k, RMSE_peach(k), RMSE_nm(k),  ...
            RMSE_golden(k), CRBth(k), MCS);
end

% ---------- 3. Gráficos de desempenho -----------------------------------
%title(sprintf('PEACH-MUSIC – %d×%d URA | L=%d | MCS=%d | SNR=%d dB', ...
 %     Mx, Mz, L, MCS, SNRdB_fixed));

% % — axes para o zoom
% ax1 = gca;
% pos = ax1.Position;
% 
% % cria eixos menores (ajuste [x y w h] ao seu gosto)
% ax2 = axes('Position',[pos(1) pos(2) 0.3 0.3]);
% 
% % plot do zoom
% semilogy(elev_vec, RMSE_peach , 'x--', ...
%          elev_vec, RMSE_nm    , 'o-', ...
%          elev_vec, RMSE_hj    , 'x--', ...
%          elev_vec, RMSE_sbplx , 'd--', ...
%          elev_vec, RMSE_golden, 's--', ...
%          elev_vec, CRBth      , 'k--','LineWidth',1.5);
% 
% xlim([0 80])
%ylim([0 0.4])

% ---------- 4. Gráficos de tempo ----------------------------------------
fig = figure('Units','centimeters','Position',[2 2 16 12]);

semilogy(elev_vec, RMSE_peach , 'x--', ...
         elev_vec, RMSE_nm    , 'o-' , ...
         elev_vec, RMSE_golden, 's--', ...
         elev_vec, CRBth      , 'k--','LineWidth',1.5);
grid on;
xlabel('SNR (dB)');
ylabel('RMSE (m)');
legend('PEACH','PEACH-NM','PEACH-Golden','CRB','Location','best');
% (opcional: título ou anotações)

% Versão com recorte máximo da figura (bounding box apertado)
exportgraphics(fig, 'peach_rmse.eps', ...
               'ContentType','vector', ...
               'BackgroundColor','none', ...
               'Resolution',600, ...
               'BoundingBox','tight');