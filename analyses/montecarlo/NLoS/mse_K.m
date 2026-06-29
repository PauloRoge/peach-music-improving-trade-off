%============================================================
%  PEACH-MUSIC – MSE em função de K_dB (NLoS aleatório)
%  Curvas: PEACH-Golden-NM com L_scatter = 1, 5, 10, 20
%  (SEM CRB: removido do parfor e do gráfico)
%============================================================
clear; clc;
startup_nlos; % define x,y,n_hiper,n_circ,deltaArea,numIterNM,tol, etc.

% ---------- Parâmetros fixos --------------------------------------------
SNRdB_ref   = 10;               % SNR fixo
K_dB_vec    = -5:1:15;           % eixo de variação de K (dB)
nK          = numel(K_dB_vec);
Lscatter_ls = [1, 5, 10, 20];   % curvas pedidas
nL          = numel(Lscatter_ls);

% ---------- Geometria do array ------------------------------------------
[URA, URA_x, URA_z, x_h, x_v, z_h, z_v] = ...
    subarrays(Mx, Mz, d_x, d_z, elev, lambda, 0);
ref = URA(1,:); % elemento de referência

% ---------- Saídas -------------------------------------------------------
MSE_golden_nm_mat = zeros(nL, nK);
time_gd_avg       = zeros(1, nK);  % tempo médio do Golden (opcional)

% ---------- Evita broadcast excessivo -----------------------------------
constL = parallel.pool.Constant(Lscatter_ls);

% ==================== Loop principal (parfor em K) =======================
parfor k = 1:nK
    KdB = K_dB_vec(k);

    mse_col   = zeros(nL,1);
    t_gd_acc  = 0;

    for iL = 1:nL
        Lscatter = constL.Value(iL);

        err2_golden_nm = 0;

        for r = 1:MCS
            % --------- Posição aleatória do usuário --------------------
            pos = [ -50 + 100*rand , 13 + 37*rand , 0 ];

            % --------- Espalhadores (posições e coeficientes) ----------
            x_sp = -50 + 100*rand(Lscatter,1);
            y_sp =   0 +  30*rand(Lscatter,1);
            z_sp = zeros(Lscatter,1);              % 2.5D: z=0
            scatterer_pos = [x_sp, y_sp, z_sp];

            modGamma   = 0.2 + 0.6*rand(Lscatter,1); % módulo 0.2–0.8
            phaseGamma = 2*pi*rand(Lscatter,1);      % fase 0–2π
            Gamma = modGamma .* exp(1j*phaseGamma);

            % ---------- Sinais NLoS -------------------------------------
            [Yh, Yv, Y] = signals_nlos_multi(pos, URA, lambda, L, alpha, ...
                                             SNRdB_ref, P_tx, Mx, Mz, KdB, ...
                                             scatterer_pos, Gamma);

            % ---------- PEACH-Golden (coarse) ---------------------------
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
            Un0 = V0(:,2:end);  % ruído (1 fonte dominante)

            est_gd_nm = nelder_mead(URA, est_gd, Un0, ...
                lambda, ref, deltaArea, numIterNM, tol, false, x, y);

            % ---------- Erro quadrático 2D ------------------------------
            err2_golden_nm = err2_golden_nm + norm(est_gd_nm(1:2) - pos(1:2))^2;
        end

        mse_col(iL) = err2_golden_nm / MCS;
    end

    MSE_golden_nm_mat(:,k) = mse_col;
    time_gd_avg(k)         = t_gd_acc / (MCS*nL);

    fprintf('K=%2d dB | MSE_G-NM: L=1:%.3f  L=5:%.3f  L=10:%.3f  L=20:%.3f  | <t>GD=%.4fs\n', ...
        KdB, mse_col(1), mse_col(2), mse_col(3), mse_col(4), time_gd_avg(k));
end

% ==================== Figura: MSE × K_dB (sem CRB) ======================
figure('Units','centimeters','Position',[2 2 16 12]); hold on;
semilogy(K_dB_vec, MSE_golden_nm_mat(1,:), 'o-' , 'LineWidth',1.5);
semilogy(K_dB_vec, MSE_golden_nm_mat(2,:), 's-' , 'LineWidth',1.5);
semilogy(K_dB_vec, MSE_golden_nm_mat(3,:), 'd-' , 'LineWidth',1.5);
semilogy(K_dB_vec, MSE_golden_nm_mat(4,:), '^-' , 'LineWidth',1.5);
hold off; grid on;
xlabel('K (dB)'); ylabel('MSE (m^2)');

% Escala log em potências de 10 (estilo do paper)
set(gca,'YScale','log');
yticks([1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3 1e4]);
yticklabels({'10^{-3}','10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}','10^{4}'});

legend({'PEACH-Golden-NM  L_{scatter}=1', ...
        'PEACH-Golden-NM  L_{scatter}=5', ...
        'PEACH-Golden-NM  L_{scatter}=10', ...
        'PEACH-Golden-NM  L_{scatter}=20'}, 'Location','best');
title(sprintf('MSE \\times K – %d\\times%d URA | L=%d | SNR=%d dB | MCS=%d', ...
    Mx, Mz, L, SNRdB_ref, MCS));

% ==================== (Opcional) Tempo × K_dB ============================
% figure('Units','centimeters','Position',[2 2 16 12]);
% semilogy(K_dB_vec, time_gd_avg, 'o-', 'LineWidth',1.5);
% grid on; xlabel('K (dB)'); ylabel('Tempo médio por realização (s)');
% legend('PEACH-Golden (tempo)','Location','best');
% title(sprintf('Tempo médio – PEACH-Golden | %d\\times%d URA | SNR=%d dB', ...
%     Mx, Mz, SNRdB_ref));
