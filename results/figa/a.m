clear; clc;
startup;
% ---------- 1. Geometria do array ---------------------------------------
[URA, URA_x, URA_z, x_h, x_v, z_h, z_v] = ...
    subarrays(Mx, Mz, d_x, d_z, elev, lambda, 0);
ref = URA(1,:);                        % elemento de referência

% ---------- 2. Varredura de SNR -----------------------------------------
SNR_dB_vec = 0:1:20;            
nsnr       = numel(SNR_dB_vec);

MSE_peach  = zeros(1,nsnr);

MSE_peach1  = zeros(1,nsnr);
MSE_peach2  = zeros(1,nsnr);
MSE_peach3  = zeros(1,nsnr);

MSE_nm     = zeros(1,nsnr);
MSE_golden = zeros(1,nsnr);
CRBth       = zeros(1,nsnr);

time_pch_avg = zeros(1,nsnr);

time_pch_avg1 = zeros(1,nsnr);
time_pch_avg2 = zeros(1,nsnr);
time_pch_avg3 = zeros(1,nsnr);

time_nm_avg  = zeros(1,nsnr);
time_gd_avg  = zeros(1,nsnr);

time_un_avg  = zeros(1,nsnr);

% ---------- 3. Loop principal (parfor em SNR) ---------------------------
parfor k = 1:nsnr
    SNRdB = SNR_dB_vec(k);

    % acumuladores de erro e tempo para este SNR
    err2_pc = 0;

    err2_pc1 = 0;
    err2_pc2 = 0;
    err2_pc3 = 0;

    err2_nm = 0;
    err2_gd = 0;
    t_pc = 0;
    t_pc1 = 0;
    t_pc2 = 0;
    t_pc3 = 0;

    t_nm = 0;
    t_gd = 0;
    t_un = 0;
    crb_sum = 0;

    for r = 1:MCS
        rng(r);                                 
        pos = [-50+100*rand,13+37*rand,0]; % pos. rand.

        % ---------- Sinais ------------------------------------------------
        [Yh, Yv, Y] = signals(pos, URA, lambda, L, ...
                              alpha, SNRdB, P_tx, Mx, Mz);

        % ---------- CRB ---------------------------------------------------
        crb_sum = crb_sum + crb(L, URA, pos, lambda, P_tx, SNRdB,2);

        % ---------- PEACH (coarse) ----------------------------------------
        tic;                                              % cronômetro
        [~, ~, est_peach] = peach( ...
            Yh, Yv, L, x, 24, x_h, z_h, x_v, z_v, ...
            ref, lambda, y, 12, pos);
        t_pc = t_pc + toc;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;                                              % cronômetro
        [~, ~, est_peach1] = peach( ...
            Yh, Yv, L, x, 96, x_h, z_h, x_v, z_v, ...
            ref, lambda, y, 48, pos);
        t_pc1 = t_pc1 + toc;
        tic;                                             
        [~, ~, est_peach2] = peach( ...
            Yh, Yv, L, x, 192, x_h, z_h, x_v, z_v, ...
            ref, lambda, y, 98, pos);
        t_pc2 = t_pc2 + toc;
        tic;                                             
        [~, ~, est_peach3] = peach( ...
            Yh, Yv, L, x, 960, x_h, z_h, x_v, z_v, ...
            ref, lambda, y, 480, pos);
        t_pc3 = t_pc3 + toc;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % ---------- Golden-PEACH (1-D scan) -----------------------------
        tic;
        [~, ~, est_golden] = peachgolden( ...
            Yh, Yv, L, x, 24, x_h, z_h, x_v, z_v, ...
            ref, lambda, y, 12, pos);
        t_gd = t_gd + toc;

        % ---------- Sub-espaço de ruído para NM e Subplex ----------------
        tic;
        Cov = (Y*Y')/L;
        [V, D] = eig(Cov); [~,idx] = sort(diag(D),'descend');
        Un = V(:, idx(2:end));                          % fonte única
        t_un = t_un + toc;

        % ---------- Nelder-Mead ------------------------------------------
        tic;
        est_nm = neldermead(URA, est_peach1, Un, lambda, ref, ...
                             deltaArea, numIterNM, tol, false, x, y);
        t_nm = t_nm + toc;

        % ---------- Acumula erros quadráticos -----------------------------
        err2_pc = err2_pc + norm(est_peach  - pos(1:2))^2;

        err2_pc1 = err2_pc1 + norm(est_peach1  - pos(1:2))^2;
        err2_pc2 = err2_pc2 + norm(est_peach2  - pos(1:2))^2;
        err2_pc3 = err2_pc3 + norm(est_peach3  - pos(1:2))^2;

        err2_nm = err2_nm + norm(est_nm     - pos(1:2))^2;
        err2_gd = err2_gd + norm(est_golden - pos(1:2))^2;
    end

    % ---------- Estatísticas por SNR -------------------------------------
    MSE_peach(k)  = err2_pc / MCS;

    MSE_peach1(k)  = err2_pc1 / MCS;
    MSE_peach2(k)  = err2_pc2 / MCS;
    MSE_peach3(k)  = err2_pc3 / MCS;
    
    MSE_nm(k)     = err2_nm / MCS;
    MSE_golden(k) = err2_gd / MCS;

    CRBth(k)       = crb_sum / MCS;

    time_un_avg(k)  = t_un / MCS;

    time_pch_avg(k) = t_pc / MCS;
    
    time_pch_avg1(k) = t_pc1 / MCS;
    time_pch_avg2(k) = t_pc2 / MCS;
    time_pch_avg3(k) = t_pc3 / MCS;

    time_nm_avg(k)  = (t_nm / MCS) + time_un_avg(k) + time_pch_avg(k);
    time_gd_avg(k)  = t_gd / MCS;

    fprintf(['SNR %2d dB | MSE→ PCH %.2f  NM %.2f' ...
             'GD %.2f | CRB %.2f  (avg de %.0f exec.)\n'], ...
            SNRdB, MSE_peach(k), MSE_nm(k), ...
            MSE_golden(k), CRBth(k), MCS);
end

%% Curvas de MSE × SNR – PEACH-MUSIC
figure('Units','centimeters','Position',[2 2 16 12]);

semilogy(SNR_dB_vec, MSE_peach , 'x--', ...
         SNR_dB_vec, MSE_peach1 , 'x--',...
         SNR_dB_vec, MSE_peach2 , 'x--',...
         SNR_dB_vec, MSE_peach3 , 'x--',...
         SNR_dB_vec, MSE_nm    , 'o-' , ...
         SNR_dB_vec, MSE_golden, 's--', ...
         SNR_dB_vec, CRBth      , 'k--','LineWidth',1.5);
grid on;
xlabel('SNR (dB)');
ylabel('MSE (m)');
legend('PEACH, N_H=24, N_C=12','PEACH, N_H=96, N_C=48','PEACH, N_H= 192, N_C=98', ...
    'PEACH, N_H=960, N_C=480','PEACH-NM, N_H=96, N_C=48','PEACH-Golden, N_H=24, N_C=12', ...
    'CRB','Location','best');
title(sprintf('PEACH-MUSIC – %d×%d URA  |  L=%d  |  MCS=%d', Mx,Mz,L,MCS));