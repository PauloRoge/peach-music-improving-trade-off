
%============================================================
% Paulo R. A. Candido Jr.
% Este script realiza análises comparativas de variantes.
%============================================================
clear; clc;
startup; % Inicializar

% -------------- 2. Modelagem Array ---------------------------
[URA, URA_x, URA_z, x_h, x_v, z_h, z_v] = subarrays(Mx, Mz, ...
    d_x, d_z, elev, lambda, 0);
ref = URA(1,:);

% -------------- 3. Loops de SNR -----------------------------
SNR_dB_vec = 0:1:20;

RMSE_peach  = zeros(size(SNR_dB_vec));
RMSE_nm     = zeros(size(SNR_dB_vec)); 
RMSE_hj     = zeros(size(SNR_dB_vec)); 
RMSE_sbplx  = zeros(size(SNR_dB_vec)); 
RMSE_golden = zeros(size(SNR_dB_vec));
CRBth       = zeros(size(SNR_dB_vec)); 
% ------------------------------------------------------------
parfor k = 1:numel(SNR_dB_vec)
    SNRdB = SNR_dB_vec(k);
    err2_pc = 0;
    err2_nm = 0;
    err2_sb = 0;
    err2_hj = 0;
    err2_gd = 0;
    crb_sum = 0;

    for r = 1:MCS
        rng(r);                        
        % ---------------------------------------------------
        % posição aleatória do usuário por realização
        x_rand = -50 + 100*rand;         % x ∈ [-50, 50]
        y_rand =  13 + 37*rand;          % y ∈ [ 13, 50]
        pos    = [x_rand, y_rand, 0];    % posição 3-D
        % ---------------------------------------------------

        % CRB desta realização
        crb_sum = crb_sum + crb(L, URA, pos, lambda, P_tx, ...
            SNRdB, alpha);

        [Yh, Yv, Y] = signals(pos, URA, lambda, L, ...
            alpha, SNRdB, P_tx, Mx, Mz);

        % ---------------------- Refinamento -------------------------
        [~, ~, est_peach] = peach_analitico(Yh, Yv, L, x, n_hiper, ...
            x_h, z_h, x_v, z_v, ref, lambda, y, n_circ, pos);

        [~, ~, est_golden] = peach_aurea( ...
        Yh, Yv, L, x, n_hiper, x_h, z_h, x_v, z_v, ref, lambda, y, n_circ, pos);
        %------------------------------------------------------------

        Un = subspace(Y, L); % Sub-espaço ruído

        % -------------------Otimização -----------------------------
        est_nelder = nelder_mead(URA, est_peach, Un, lambda, ...
            ref, deltaArea, numIterNM, 1e-6, true, x, y);

        hj_est = hooke_jeeves(URA, est_peach, Un, lambda, ...
            ref, deltaArea, 100, tol, x, y, true);

        est_sbplx = subplex_wrapper(URA, est_peach, Un, lambda, ref, x, y, tol, 50)
        % -----------------------------------------------------------
        
        % ------------- Erro quadrático acumulado -------------------
        err2_pc  = err2_pc  + norm(est_peach  - pos(1:2))^2;
        err2_nm     = err2_nm     + norm(est_nelder - pos(1:2))^2;
        err2_hj     = err2_hj     + norm(hj_est     - pos(1:2))^2;
        err2_sb  = err2_sb  + norm(est_sbplx  - pos(1:2))^2;
        err2_gd = err2_gd + norm(est_golden - pos(1:2))^2;
        %------------------------------------------------------------
    end

    % ----------------- RMSE médio para este SNR --------------------
    RMSE_peach(k)  = sqrt(err2_pc / MCS);
    RMSE_nm(k)     = sqrt(err2_nm / MCS);
    RMSE_hj(k)     = sqrt(err2_hj / MCS);
    RMSE_sbplx(k)  = sqrt(err2_sb / MCS);
    RMSE_golden(k) = sqrt(err2_gd / MCS);
    %----------------------------------------------------------------

    % ----------------------- CRB médio -----------------------------
    CRBth(k) = crb_sum / MCS;

    fprintf(['SNR = %2d dB | PCH=%.2f m  NM=%.2f m ' ...
        ' HJ =%.2f m  sbplx =%.2f m GD=%.3f  CRB = %.3f m \n'], ...
        SNRdB, RMSE_peach(k), RMSE_nm(k), RMSE_hj(k), ...
        RMSE_sbplx(k),RMSE_golden(k), CRBth(k));
end

% ---------------------- 4. PLOT ------------------------------------
figure;
semilogy(SNR_dB_vec, RMSE_peach,  'x--', 'LineWidth', 1.5); hold on;
semilogy(SNR_dB_vec, RMSE_nm,     'o-',  'LineWidth', 1.5);
semilogy(SNR_dB_vec, RMSE_hj,     'x--', 'LineWidth', 1.5);
semilogy(SNR_dB_vec, RMSE_sbplx,  'd--', 'LineWidth', 1.5);
semilogy(SNR_dB_vec, RMSE_golden, 's--', 'LineWidth', 1.5);
semilogy(SNR_dB_vec, CRBth,       'k--', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Erro (m)');
legend('PEACH','PEACH + Nelder Mead','PEACH + Hooke Jeeves', ...
       'PEACH + Subplex','Golden PEACH','CRB');
title(sprintf('PEACH-MUSIC - %d×%d URA - L=%d - MCS=%d', ...
       Mx, Mz, L, MCS));

% ---------------- FUNÇÃO AUXILIAR ----------------------------
function Un = subspace(Y, L)
    Cov = (Y * Y') / L;
    [eigenvectors,eigenvalues] = eig(Cov);
    [~, idx] = sort(diag(eigenvalues),'descend');
    eigenvectors = eigenvectors(:, idx);
    Un = eigenvectors(:, 2:end);
end
