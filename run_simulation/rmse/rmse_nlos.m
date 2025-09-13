%============================================================
%  PEACH-MUSIC – Comparação de variantes com PARFOR
%============================================================
clear; clc;
startup_nlos; % parâmetros gerais do projeto

K_dB = 3;

scatterer_pos = [  % Lx3
    -35, 15, 2;
    -8,  10, 1.2;
    20, 5,  2.5
];

Gamma = [ ...
    0.5*exp(1j*pi/3);
    0.3*exp(1j*1.1*pi);
    0.4*exp(1j*0.2*pi)
];
% ---------- 1. Geometria do array ---------------------------------------
[URA, URA_x, URA_z, x_h, x_v, z_h, z_v] = ...
    subarrays(Mx, Mz, d_x, d_z, elev, lambda, 0);
ref = URA(1,:); % elemento de referência

% ---------- 2. Varredura de SNR -----------------------------------------
SNR_dB_vec = 0:1:20;                   % passo de 2 dB para reduzir tempo
nsnr       = numel(SNR_dB_vec);

%% ---------------- Pré-alocações -----------------------------------------
MSE_golden_los  = zeros(1,nsnr);
MSE_golden_NM_los = zeros(1,nsnr);
MSE_golden_nlos = zeros(1,nsnr);
MSE_golden_NM_nlos = zeros(1,nsnr);
CRBth = zeros(1,nsnr);

time_gd_los_avg  = zeros(1,nsnr);
time_gd_nlos_avg = zeros(1,nsnr);

%% ---------------- Loop principal (parfor em SNR) ------------------------
parfor k = 1:nsnr
    SNRdB = SNR_dB_vec(k);

    % Acumuladores por SNR
    err2_los = 0;
    err2_los_nm = 0;
    err2_nlos_nm = 0;
    err2_nlo = 0;
    t_los    = 0;
    t_nlo    = 0;
    crb_sum  = 0;

    for r = 1:MCS
        % Semente reprodutível (única por (k,r))
        rng(10^6*k + r,'twister');

        % Posição aleatória do usuário (ajuste conforme seu cenário)
        pos = [ -50 + 100*rand , 13 + 37*rand , 0 ];

        % ---------- Sinais LoS ------------------------------------------
        % Saída esperada: Wh, Wv, Y_los (ou W) compatível com seu peach_golden
        [Wh, Wv, W] = signals_los(pos, URA, lambda, L, alpha, SNRdB, P_tx, Mx, Mz);

        % ---------- Sinais NLoS -----------------------------------------
        % Saída esperada: Yh, Yv, Y_nlos
        [Yh, Yv, Y] = signals_nlos_multi(pos, URA, lambda, L, alpha, SNRdB, ...
                                              P_tx, Mx, Mz, K_dB, scatterer_pos, Gamma);

        % ---------- CRB (único usuário) ---------------------------------
        crb_sum = crb_sum + crb(L, URA, pos, lambda, P_tx, SNRdB, 2);

        % ---------- PEACH-Golden (LoS) ----------------------------------
        t0 = tic;
        [~, ~, est_gd_los] = peach_golden(Wh, Wv, L, x, n_hiper, x_h, z_h, x_v, z_v, ...
                                          ref, lambda, y, n_circ, pos);
        t_los = t_los + toc(t0);

        % ---------- PEACH-Golden (NLoS) ---------------------------------
        t0 = tic;
        [~, ~, est_gd_nlo] = peach_golden(Yh, Yv, L, x, n_hiper, x_h, z_h, x_v, z_v, ...
                                          ref, lambda, y, n_circ, pos);
        t_nlo = t_nlo + toc(t0);
        
        Cov_nlos = (Y*Y')/L;
        [V,D] = eig(Cov_nlos); [~,idx] = sort(diag(D),'descend');
        Un_nlos = V(:, idx(2:end));
        
        Cov_los = (W*W')/L;
        [V,D] = eig(Cov_los); [~,idx] = sort(diag(D),'descend');
        Un_los = V(:, idx(2:end));
        
        est_nm_los = nelder_mead(URA, est_gd_los, Un_los, lambda, ref, ...
            deltaArea, numIterNM, tol, false, x, y);
        
        est_nm_nlos = nelder_mead(URA, est_gd_nlo, Un_nlos, lambda, ref, ...
            deltaArea, numIterNM, tol, false, x, y);

        % ---------- Erros quadráticos 2D (x,y) --------------------------
        err2_los = err2_los + norm(est_gd_los(1:2) - pos(1:2))^2;
        err2_nlo = err2_nlo + norm(est_gd_nlo(1:2) - pos(1:2))^2;
        err2_los_nm = err2_los_nm + norm(est_nm_los(1:2) - pos(1:2))^2;
        err2_nlos_nm = err2_nlos_nm + norm(est_nm_nlos(1:2) - pos(1:2))^2;
    end

    % Médias por SNR
    MSE_golden_los(k)  = err2_los / MCS;
    MSE_golden_nlos(k) = err2_nlo / MCS;
    MSE_golden_nm_los(k)  = err2_los_nm / MCS;
    MSE_golden_nm_nlos(k) = err2_nlos_nm / MCS;
    CRBth(k)           = crb_sum  / MCS;

    time_gd_los_avg(k)  = t_los / MCS;
    time_gd_nlos_avg(k) = t_nlo / MCS;

    fprintf('SNR %2d dB | MSE  LoS %.3f | MSE NLoS %.3f | CRB %.3f | <t> LoS %.4fs | <t> NLoS %.4fs\n', ...
        SNRdB, MSE_golden_los(k), MSE_golden_nlos(k), CRBth(k), ...
        time_gd_los_avg(k), time_gd_nlos_avg(k));
end

%% ----------------- Curvas de MSE × SNR ----------------------------------
figure('Units','centimeters','Position',[2 2 16 12]);
semilogy(SNR_dB_vec, MSE_golden_los , 'o-', ...
         SNR_dB_vec, MSE_golden_nlos , 'o-', ...
         SNR_dB_vec, MSE_golden_nm_los, 's-', ...
         SNR_dB_vec, MSE_golden_nm_nlos, 's-', ...
         SNR_dB_vec, CRBth          , 'k--', 'LineWidth',1.5);
grid on; xlabel('SNR (dB)'); ylabel('MSE (m^2)');
legend('Golden (LoS)','Golden (NLoS)','Golden-NM (LoS)','Golden-NM (NLoS)','CRB','Location','best');
title(sprintf('PEACH-Golden – %d\\times%d URA  |  L=%d  |  MCS=%d', Mx, Mz, L, MCS));

%% ----------------- Curvas de Tempo × SNR --------------------------------
figure('Units','centimeters','Position',[2 2 16 12]);
semilogy(SNR_dB_vec, time_gd_los_avg , 'o-', ...
         SNR_dB_vec, time_gd_nlos_avg, 's-', 'LineWidth',1.5);
grid on; xlabel('SNR (dB)'); ylabel('Tempo médio por realização (s)');
legend('Golden (LoS)','Golden (NLoS)','Location','best');
title(sprintf('Tempo médio – PEACH-Golden  |  %d\\times%d URA', Mx, Mz));
