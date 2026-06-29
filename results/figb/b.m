clear; clc;
startup;
% ---------- 1. Configurações --------------------------------------------
SNRdB_fixed = 15;            % SNR constante (dB)
elev_vec    = 0:5:100;       % elevação do array (m) — passo de 5 m
nelev       = numel(elev_vec);

RMSE_peach  = zeros(1,nelev);
RMSE_nm     = zeros(1,nelev);
RMSE_golden = zeros(1,nelev);
CRBth       = zeros(1,nelev);

time_pch_avg = zeros(1,nelev);
time_nm_avg  = zeros(1,nelev);
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
    t_gd = 0; t_un = 0;
    crb_sum = 0;

    for r = 1:MCS
        rng(r);                                   
        pos = [ -50 + 100*rand , 13 + 37*rand , 0 ];

        % ---------- Sinais ----------------------------------------------
        [Yh, Yv, Y] = signals(pos, URA, lambda, L, ...
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
        [~, ~, est_golden] = peachgolden( ...
            Yh, Yv, L, x, n_hiper, x_h, z_h, x_v, z_v, ...
            ref, lambda, y, n_circ, pos);
        t_gd = t_gd + toc;

        % ---------- Sub-espaço de ruído ---------------------------------
        tic;
        Cov = (Y*Y')/L;
        [V,D] = eig(Cov); [~,idx] = sort(diag(D),'descend');
        Un = V(:, idx(2:end));                       % fonte única
        t_un = t_un + toc;

        % ---------- Nelder-Mead -----------------------------------------
        tic;
        est_nm = neldermead(URA, est_peach, Un, lambda, ref, ...
                             deltaArea, numIterNM, tol, false, x, y);
        t_nm = t_nm + toc;

        % ---------- Acumula erros ---------------------------------------
        err2_pc = err2_pc + norm(est_peach  - pos(1:2))^2;
        err2_nm = err2_nm + norm(est_nm     - pos(1:2))^2;
        err2_gd = err2_gd + norm(est_golden - pos(1:2))^2;
    end

    % ---------- Estatísticas por elevação -------------------------------
    RMSE_peach(k)  = err2_pc / MCS;
    RMSE_nm(k)     = err2_nm / MCS;
    RMSE_golden(k) = err2_gd / MCS;
    CRBth(k)       = crb_sum / MCS;

    time_un_avg(k)  = t_un / MCS;
    time_pch_avg(k) = t_pc / MCS;
    time_nm_avg(k)  = (t_nm / MCS) + time_un_avg(k) + time_pch_avg(k);
    time_gd_avg(k)  = t_gd / MCS;

    fprintf(['Elev %3.0f m | RMSE→ PCH %.2f  NM %.2f  HJ %.2f  SB %.2f  ' ...
             'GD %.2f | CRB %.2f  (avg de %.0f exec.)\n'], ...
            elev_k, RMSE_peach(k), RMSE_nm(k),  ...
            RMSE_golden(k), CRBth(k), MCS);
end