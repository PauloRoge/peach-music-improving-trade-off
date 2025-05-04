clc;
startup;

% ------------ VARREDURA EM ELEVAÇÃO -------------------
elev_vec = 0:1:100;  
RMSE_nm    = zeros(size(elev_vec));  % RMSE do Nelder-Mead
RMSE_peach = zeros(size(elev_vec));  % RMSE do PEACH puro
CRBth      = zeros(size(elev_vec));

for i = 1:numel(elev_vec)
    elev = elev_vec(i);

    [URA, ~, ~, x_h, x_v, z_h, z_v] = subarrays(Mx, Mz, d_x, d_z, elev, lambda, 0);
    ref = URA(1,:);
    CRBth(i) = crb(L, URA, pos, lambda, P_tx, SNR_dB, alpha);

    err2_nm    = 0;
    err2_peach = 0;

    for r = 1:MCS
    
        [Yh, Yv, Y] = signals(pos, URA, lambda, L, alpha, SNR_dB, P_tx, Mx, Mz);
        %[Yh, Yv, Y] = signals_fig1(UEs, URA, lambda, L, SNR_dB, P_tx, Mx, Mz);

        
        % Estimativa PEACH puro
        [~, ~, estPos] = peach_analitico(Yh, Yv, L, x, n_hiper, ...
            x_h, z_h, x_v, z_v, ref, lambda, y, n_circ, pos);
        err2_peach = err2_peach + norm(estPos - pos(1:2))^2;

        % Refinamento Nelder-Mead
        Cov = (Y * Y') / L;
        [V, D] = eig(Cov);
        [~, j] = sort(diag(D), 'descend');
        Un = V(:, j(2:end));

        % [nm_est, ~] = nelder_mead(URA, estPos, Un, lambda, ref, ...
        %     deltaArea, numIterNM, 1e-6, true, x, y);

        [subplex_est, hist] = subplex_wrapper(URA, estPos, Un, lambda, ref, x, y, 1e-5, 20);

        err2_nm = err2_nm + norm(subplex_est - pos(1:2))^2;
    end

    RMSE_peach(i) = sqrt(err2_peach / MCS);
    RMSE_nm(i)    = sqrt(err2_nm    / MCS);

    fprintf('Elev = %3d m | RMSE_peach = %.3f m | RMSE_NM = %.3f m | CRB = %.3f m\n', ...
            elev, RMSE_peach(i), RMSE_nm(i), CRBth(i));
end

% ---------------- PLOT RESULTADO -------------------------
figure;
semilogy(elev_vec, RMSE_peach, 'b-', 'LineWidth', 2); hold on;
semilogy(elev_vec, RMSE_nm,    '-',   'LineWidth', 2);
semilogy(elev_vec, CRBth,      'k--', 'LineWidth', 1.5);
xlabel('Elevação do array (m)');
ylabel('Erro (m)');
legend('PEACH puro','PEACH + Subplex','CRB teórico','Location','southwest');
grid on;
title(sprintf('PEACH-MUSIC  |  SNR = %d dB,  %d×%d URA,  L = %d,  MCS = %d', ...
       SNR_dB, Mx, Mz, L, MCS));
