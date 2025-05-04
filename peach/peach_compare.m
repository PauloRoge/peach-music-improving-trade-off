%============================================================
% Paulo R. A. Candido Jr.
% Este script realiza análises comparativas de variantes.
%============================================================
clear; clc;

% INICIALIZAR
startup;

% -------------- 2. PREPARA ARRAY ----------------------------
[URA, URA_x, URA_z, x_h, x_v, z_h, z_v] = subarrays(Mx, Mz, d_x, d_z, elev, lambda, 0);
ref = URA(1,:);

% -------------- 3. LAÇO DE SNRs -----------------------------
SNR_dB_vec      = 0:1:20;   % vetores de SNR

RMSE_peach      = zeros(size(SNR_dB_vec));  % PEACH_analítico
RMSE_peach2 = zeros(size(SNR_dB_vec));  % PEACH_2
CRBth           = zeros(size(SNR_dB_vec));

for k = 1:numel(SNR_dB_vec)
    SNRdB = SNR_dB_vec(k);
    err2_analitico = 0;
    err2_test      = 0;

    % CRB para este SNR
    CRBth(k) = crb(L, URA, pos, lambda, P_tx, SNRdB, alpha);

    for r = 1:MCS
        % Gera nova posição aleatória para cada realização
        %rng(r);
        % x_rand = x(1) + (x(2) - x(1)) * rand;
        % y_rand = y(1) + (y(2) - y(1)) * rand;
        % pos = [x_rand, y_rand, 0];

        [Yh, Yv, Y] = signals(pos, URA, lambda, L, alpha, SNRdB, P_tx, Mx, Mz);

        % PEACH analítico
        [~, ~, est_peach] = peach_analitico(Yh, Yv, L, x, n_hiper, ...
            x_h, z_h, x_v, z_v, ref, lambda, y, n_circ, pos);
        
        %PEACH (TESTE) ====================================================
        [~, ~, est_aleatorio] = randomized_peach(Yh, Yv, L, x, n_hiper, ...
            x_h, z_h, x_v, z_v, ref, ...
            lambda, y, n_circ, pos);
        % =================================================================

        [Un] = subspace(Y, L); % divisao do subespaco
        %===============================================================
        %                           OTIMIZAÇÃO
        % ==============================================================
        % % Otimização peach referencia
        % [est_nelder, ~] = nelder_mead(URA, est_peach, Un, lambda, ...
        %      ref, deltaArea, numIterNM, 1e-6, true, x, y);
        % 
        % % Otimização peach teste 
        % [subplex_est, hist] = subplex_wrapper(URA, est_peach, Un, ...
        %     lambda, ref, x, y, 1e-6, 200);
        % % ============================================================

        err2_analitico = err2_analitico + norm(est_peach - pos(1:2))^2;
        err2_test = err2_test + norm(est_aleatorio  - pos(1:2))^2;
    end

    RMSE_peach(k) = sqrt(err2_analitico / MCS);
    RMSE_peach2(k) = sqrt(err2_test / MCS);

    % fprintf(['SNR = %3d dB | RMSE_peach = %.3f m |' ...
    %     ' RMSE_teste = %.3f m\n '], ...
    %         SNRdB, RMSE_peach(k), RMSE_peach2(k));

    fprintf(['SNR = %+3d dB | RMSE_peach = %.3f m |' ...
        ' RMSE_test = %.3f m | CRB = %.3f m\n'], ...
        SNRdB, RMSE_peach(k), RMSE_peach2(k), CRBth(k));
end

% -------------- 4. PLOT ----------------------------------------
figure; 
semilogy(SNR_dB_vec, RMSE_peach, 'o-',  'LineWidth', 1.5); hold on;
semilogy(SNR_dB_vec, RMSE_peach2, 'x--', 'LineWidth', 1.5);
semilogy(SNR_dB_vec, CRBth, 's--', 'LineWidth', 1.5);

grid on;
xlabel('SNR (dB)');
ylabel('Erro (m)');
legend('PEACH analítico','PEACH test','CRB teórico','Location','southwest');
title(sprintf('PEACH-MUSIC  |  %d×%d URA,  L = %d,  N_{real} = %d', ...
       Mx, Mz, L, MCS));

function [Un] = subspace(Y, L)
        Cov = (Y * Y') / L;
        [eigenvectors, eigenvalues] = eig(Cov); 
        estimated_sources = 1;
        [~, i] = sort(diag(eigenvalues), 'descend'); 
        eigenvectors = eigenvectors(:, i);
        Un = eigenvectors(:, estimated_sources+1:end);
end
