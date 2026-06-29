%============================================================
%  PEACH-MUSIC – Comparação (NLoS aleatório) com PARFOR
%  L_scatter = 10 (fixo)
%  K_dB ∈ {-5, 5, 15}
%  Curvas: PEACH, PEACH-Golden, PEACH-Golden+NM  +  CRB_nlos
%============================================================
clear; clc;
startup_nlos;  % define x,y,n_hiper,n_circ,deltaArea,numIterNM,tol, etc.

% ---------- Parâmetros do cenário ---------------------------------------
Lscatter = 10;                         % nº de espalhadores (fixo)
Kset_dB  = [-5, 5, 10];                % valores de K (dB)
nK       = numel(Kset_dB);

% ---------- Geometria do array ------------------------------------------
[URA, URA_x, URA_z, x_h, x_v, z_h, z_v] = ...
    subarrays(Mx, Mz, d_x, d_z, elev, lambda, 0);
ref = URA(1,:);

% ---------- Eixo de SNR --------------------------------------------------
SNR_dB_vec = 0:1:20;
nsnr       = numel(SNR_dB_vec);

% ---------- Pré-alocações ------------------------------------------------
MSE_peach_mat      = zeros(nK, nsnr);
MSE_golden_mat     = zeros(nK, nsnr);
MSE_golden_nm_mat  = zeros(nK, nsnr);
CRB_mat            = zeros(nK, nsnr);    % CRB_nlos por K
time_gd_avg_mat    = zeros(nK, nsnr);    % (opcional) tempo do Golden

% ---------- Evita broadcast do vetor de K --------------------------------
constK = parallel.pool.Constant(Kset_dB);

% ====================== Loop principal (parfor em SNR) ===================
parfor ksnr = 1:nsnr
    SNRdB = SNR_dB_vec(ksnr);

    mse_peach_col     = zeros(nK,1);
    mse_golden_col    = zeros(nK,1);
    mse_golden_nm_col = zeros(nK,1);
    crb_col           = zeros(nK,1);
    tgd_col           = zeros(nK,1);

    Klist = constK.Value;  % cópia local no worker

    for ik = 1:nK
        K_dB = Klist(ik);

        err2_peach      = 0;
        err2_golden     = 0;
        err2_golden_nm  = 0;
        crb_sum         = 0;
        t_gd_acc        = 0;

        for r = 1:MCS
            % --------- Usuário aleatório --------------------------------
            pos = [ -50 + 100*rand , 13 + 37*rand , 0 ];

            % --------- Espalhadores (posições e coeficientes) -----------
            x_sp = -50 + 100*rand(Lscatter,1);
            y_sp =   0 +  30*rand(Lscatter,1);
            z_sp = zeros(Lscatter,1);                 % 2.5D: z=0
            scatterer_pos = [x_sp, y_sp, z_sp];

            modGamma   = 0.2 + 0.6*rand(Lscatter,1);  % módulo 0.2–0.8
            phaseGamma = 2*pi*rand(Lscatter,1);       % fase 0–2π
            Gamma = modGamma .* exp(1j*phaseGamma);

            % ---------- Sinais NLoS -------------------------------------
            [Yh, Yv, Y] = signals_nlos_multi(pos, URA, lambda, L, alpha, ...
                                             SNRdB, P_tx, Mx, Mz, K_dB, ...
                                             scatterer_pos, Gamma);

            % ---------- PEACH (baseline) --------------------------------
            [~, ~, est_peach] = peach( ...
                Yh, Yv, L, x, 192, x_h, z_h, x_v, z_v, ...
                ref, lambda, y, 96, pos);

            % ---------- PEACH-Golden (coarse com GSS) -------------------
            t0 = tic;
            [~, ~, est_gd] = peach_golden( ...
                Yh, Yv, L, x, n_hiper, x_h, z_h, x_v, z_v, ...
                ref, lambda, y, n_circ, pos);
            t_gd_acc = t_gd_acc + toc(t0);

            % ---------- NM (refino) -------------------------------------
            Cov0 = (Y*Y')/L;
            [V0,D0] = eig(Cov0);
            [~,ord0] = sort(diag(D0),'descend');
            V0  = V0(:,ord0);
            Un0 = V0(:,2:end);                         % ruído (1 fonte)

            est_gd_nm = nelder_mead(URA, est_gd, Un0, ...
                lambda, ref, deltaArea, numIterNM, tol, false, x, y);

            % ---------- Erros quadráticos 2D ----------------------------
            err2_peach      = err2_peach      + norm(est_peach(1:2) - pos(1:2))^2;
            err2_golden     = err2_golden     + norm(est_gd(1:2)    - pos(1:2))^2;
            err2_golden_nm  = err2_golden_nm  + norm(est_gd_nm(1:2) - pos(1:2))^2;

            % ---------- CRB NLoS (2D) -----------------------------------
            crb_sum = crb_sum + crb_nlos(L, URA, pos, lambda, P_tx, ...
                                         SNRdB, 2, K_dB, scatterer_pos, Gamma);
        end

        mse_peach_col(ik)     = err2_peach     / MCS;
        mse_golden_col(ik)    = err2_golden    / MCS;
        mse_golden_nm_col(ik) = err2_golden_nm / MCS;
        crb_col(ik)           = crb_sum        / MCS;
        tgd_col(ik)           = t_gd_acc       / MCS;
    end

    MSE_peach_mat(:,ksnr)     = mse_peach_col;
    MSE_golden_mat(:,ksnr)    = mse_golden_col;
    MSE_golden_nm_mat(:,ksnr) = mse_golden_nm_col;
    CRB_mat(:,ksnr)           = crb_col;
    time_gd_avg_mat(:,ksnr)   = tgd_col;

    fprintf(['SNR %2d dB | MSE[PEACH/G/G+NM] (K=-5,5,15): ' ...
             '[%.3f %.3f %.3f] / [%.3f %.3f %.3f] / [%.3f %.3f %.3f] | ' ...
             'CRB: [%.3f %.3f %.3f]\n'], ...
        SNRdB, ...
        mse_peach_col(1),  mse_peach_col(2),  mse_peach_col(3), ...
        mse_golden_col(1), mse_golden_col(2), mse_golden_col(3), ...
        mse_golden_nm_col(1), mse_golden_nm_col(2), mse_golden_nm_col(3), ...
        crb_col(1), crb_col(2), crb_col(3));
end

% ====================== Plot: MSE × SNR (9 + 3 curvas) ===================
figure('Units','centimeters','Position',[2 2 16 12]); hold on;

% ---- PEACH
semilogy(SNR_dB_vec, MSE_peach_mat(1,:), 'o-' , 'LineWidth',1.5);
semilogy(SNR_dB_vec, MSE_peach_mat(2,:), 'o--', 'LineWidth',1.5);
semilogy(SNR_dB_vec, MSE_peach_mat(3,:), 'o:' , 'LineWidth',1.5);

% ---- PEACH-Golden
semilogy(SNR_dB_vec, MSE_golden_mat(1,:), 's-' , 'LineWidth',1.5);
semilogy(SNR_dB_vec, MSE_golden_mat(2,:), 's--', 'LineWidth',1.5);
semilogy(SNR_dB_vec, MSE_golden_mat(3,:), 's:' , 'LineWidth',1.5);

% ---- PEACH-Golden + NM
semilogy(SNR_dB_vec, MSE_golden_nm_mat(1,:), 'd-' , 'LineWidth',1.5);
semilogy(SNR_dB_vec, MSE_golden_nm_mat(2,:), 'd--', 'LineWidth',1.5);
semilogy(SNR_dB_vec, MSE_golden_nm_mat(3,:), 'd:' , 'LineWidth',1.5);

% ---- CRB_nlos (linhas pretas tracejadas)
semilogy(SNR_dB_vec, CRB_mat(1,:), 'k-.' , 'LineWidth',1.2);
semilogy(SNR_dB_vec, CRB_mat(2,:), 'k--' , 'LineWidth',1.2);
semilogy(SNR_dB_vec, CRB_mat(3,:), 'k:'  , 'LineWidth',1.4);

hold off; grid on;
xlabel('SNR (dB)'); ylabel('MSE / CRB (m^2)');

% Escala log em potências de 10
set(gca,'YScale','log');
yticks([1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3 1e4]);
yticklabels({'10^{-3}','10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}','10^{4}'});

lgd = legend({ ...
 'PEACH K=-5','PEACH K=5','PEACH K=15', ...
 'PEACH-Golden K=-5','PEACH-Golden K=5','PEACH-Golden K=15', ...
 'PEACH-Golden+NM K=-5','PEACH-Golden+NM K=5','PEACH-Golden+NM K=15', ...
 'CRB NLoS K=-5','CRB NLoS K=5','CRB NLoS K=15'}, ...
 'Location','north','NumColumns',2);

lgd.FontSize = 8;   % ajuste conforme necessário

title(sprintf('MSE & CRB \\times SNR – %d\\times%d URA | L=%d | L_{scatter}=10 | MCS=%d', ...
    Mx, Mz, L, MCS));

% ---------------------- (Opcional) Tempo × SNR ---------------------------
% figure('Units','centimeters','Position',[2 2 16 12]); hold on;
% semilogy(SNR_dB_vec, time_gd_avg_mat(1,:), 'o-', 'LineWidth',1.5);
% semilogy(SNR_dB_vec, time_gd_avg_mat(2,:), 'o--','LineWidth',1.5);
% semilogy(SNR_dB_vec, time_gd_avg_mat(3,:), 'o:','LineWidth',1.5);
% hold off; grid on; xlabel('SNR (dB)'); ylabel('<tempo GSS> (s)');
% legend('K=-5','K=5','K=15','Location','best');
% title('Tempo médio – PEACH-Golden');
