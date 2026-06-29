% peach_op_comp_time.m
clear; clc;
startup;   % carrega parâmetros: freq, lambda, L, P_tx, alpha, MCS, SNR_dB, Mx, Mz, d_x, d_z, elev, x_grid, y_grid, n_hiper, n_circ, numIterNM, tol, deltaArea

%% 1. RMSE vs SNR (para elev fixed de startup)
[URA, URA_x, URA_z, x_h, x_v, z_h, z_v] = subarrays(Mx, Mz, d_x, d_z, elev, lambda, false);
ref = URA(1,:);
SNR_dB_vec = 0:1:20;
RMSE_peach  = zeros(size(SNR_dB_vec));
RMSE_nm     = zeros(size(SNR_dB_vec)); 
RMSE_hj     = zeros(size(SNR_dB_vec)); 
RMSE_sbplx  = zeros(size(SNR_dB_vec)); 
RMSE_golden = zeros(size(SNR_dB_vec));
CRBth       = zeros(size(SNR_dB_vec));

parfor k = 1:numel(SNR_dB_vec)
    SNRdB = SNR_dB_vec(k);
    err2_pc = 0; err2_nm = 0; err2_hj = 0; err2_sb = 0; err2_gd = 0;
    crb_sum = 0;
    for r = 1:MCS
        rng(r);
        x_rand = -50 + 100*rand;
        y_rand =  13 + 37*rand;
        pos = [x_rand, y_rand, 0];
        crb_sum = crb_sum + crb(L, URA, pos, lambda, P_tx, SNRdB,2);
        [Yh, Yv, Y] = signals(pos, URA, lambda, L, alpha, SNRdB, P_tx, Mx, Mz);

        [~,~,est_pc] = peach_analitico(Yh, Yv, L, x_grid, n_hiper, x_h, z_h, x_v, z_v, ref, ...
            lambda, y_grid, n_circ, pos);

        [~,~,est_gd] = peach_aurea(Yh, Yv, L, x_grid, n_hiper, x_h, z_h, x_v, z_v, ref, ...
            lambda, y_grid, n_circ, pos);

        Un = subspace(Y, L);

        hj_est = hooke_jeeves(URA, est_pc, Un, lambda, ref, deltaArea, 15, tol, x_grid, y_grid, true);

        [est_nm,~] = nelder_mead(URA, est_pc, Un, lambda, ref, deltaArea, numIterNM, ...
            tol, false, x_grid, y_grid);

        [est_sb,~]   = subplex_wrapper(URA, est_pc, Un, lambda, ...
                                      ref, x, y, tol, numIterNM);

        err2_pc = err2_pc + norm(est_pc - pos(1:2))^2;
        err2_gd = err2_gd + norm(est_gd - pos(1:2))^2;
        err2_hj = err2_hj + norm(hj_est - pos(1:2))^2;
        err2_nm = err2_nm + norm(est_nm - pos(1:2))^2;
        err2_sb = err2_sb + norm(est_sb - pos(1:2))^2;
    end
    RMSE_peach(k)  = sqrt(err2_pc / MCS);
    RMSE_golden(k) = sqrt(err2_gd / MCS);
    RMSE_hj(k)     = sqrt(err2_hj / MCS);
    RMSE_nm(k)     = sqrt(err2_nm / MCS);
    RMSE_sbplx(k)  = sqrt(err2_sb / MCS);
    CRBth(k)       = crb_sum / MCS;
end

% Plot RMSE vs SNR
figure;
semilogy(SNR_dB_vec, RMSE_peach,  'x--','LineWidth',1.5); hold on;
semilogy(SNR_dB_vec, RMSE_nm,     'o-','LineWidth',1.5);
semilogy(SNR_dB_vec, RMSE_hj,     's--','LineWidth',1.5);
semilogy(SNR_dB_vec, RMSE_sbplx,  'd--','LineWidth',1.5);
semilogy(SNR_dB_vec, RMSE_golden, 'v--','LineWidth',1.5);
semilogy(SNR_dB_vec, CRBth,       'k--','LineWidth',1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Erro de Posição (m)');
legend('PEACH','PEACH+NM','PEACH+HJ','PEACH+Subplex','Golden','CRB','Location','best');
title(sprintf('PEACH-MUSIC – %dx%d URA – elev = %d m – MCS = %d', Mx, Mz, elev, MCS));

%% 2. Tempo médio de execução vs elevação (para SNR fixo de startup)
elev_vec   = 0:10:50; 
time_pch   = zeros(size(elev_vec));
time_nm    = zeros(size(elev_vec));
time_hj    = zeros(size(elev_vec));
time_sb    = zeros(size(elev_vec));
time_gd    = zeros(size(elev_vec));

for ei = 1:numel(elev_vec)
    elev_i = elev_vec(ei);
    [URA, URA_x, URA_z, x_h, x_v, z_h, z_v] = subarrays(Mx, Mz, d_x, d_z, elev_i, lambda, false);
    ref = URA(1,:);
    t_p=0; t_nm=0; t_hj=0; t_sb=0; t_gd=0;
    for r = 1:MCS
        rng(r);
        x_rand = -50 + 100*rand;
        y_rand =  13 + 37*rand;
        pos = [x_rand, y_rand, 0];
        [Yh, Yv, Y] = signals(pos, URA, lambda, L, alpha, SNR_dB, P_tx, Mx, Mz);

        tic; peach_analitico(Yh, Yv, L, x_grid, n_hiper, x_h, z_h, x_v, z_v, ref, lambda, y_grid, n_circ, pos);
        t_p = t_p + toc;

        tic; peach_aurea(Yh, Yv, L, x_grid, n_hiper, x_h, z_h, x_v, z_v, ref, lambda, y_grid, n_circ, pos);
        t_gd = t_gd + toc;

        Un = subspace(Y, L);

        tic; hooke_jeeves(URA, pos, Un, lambda, ref, deltaArea, 15, tol, x_grid, y_grid, false);
        t_hj = t_hj + toc;

        tic; nelder_mead(URA, pos, Un, lambda, ref, deltaArea, numIterNM, tol, false, x_grid, y_grid);
        t_nm = t_nm + toc;

        tic; subplex_wrapper(URA, pos, Un, lambda, ...
                                      ref, x, y, tol, numIterNM);
        t_sb = t_sb + toc;
    end
    time_pch(ei) = t_p / MCS;
    time_gd (ei) = t_gd / MCS;
    time_hj(ei)   = t_hj / MCS;
    time_nm(ei)   = t_nm / MCS;
    time_sb(ei)   = t_sb / MCS;
end

% Plot tempo vs elevação
figure;
semilogy(elev_vec, time_pch, 'o-','LineWidth',2); hold on;
semilogy(elev_vec, time_nm,  's-','LineWidth',2);
semilogy(elev_vec, time_hj,  'd-','LineWidth',2);
semilogy(elev_vec, time_sb,  'v-','LineWidth',2);
semilogy(elev_vec, time_gd,  'x-','LineWidth',2);
grid on;
xlabel('Elevação do array (m)');
ylabel('Tempo médio de execução (s)');
legend('PEACH','Nelder-Mead','Hooke-Jeeves','Subplex','Golden','Location','best');
title(sprintf('Tempo de Execução | SNR = %d dB – %dx%d URA – L = %d – MCS = %d', ...
      SNR_dB, Mx, Mz, L, MCS));

%% Função auxiliar (já vinha no peach_op_comp original)
function Un = subspace(Y, L)
    Cov = (Y * Y') / L;
    [V, D] = eig(Cov);
    [~, idx] = sort(diag(D),'descend');
    V = V(:, idx);
    Un = V(:, 2:end);
end
