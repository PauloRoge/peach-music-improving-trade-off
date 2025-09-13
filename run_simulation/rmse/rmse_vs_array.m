% ==============================================================
%  –– compara PEACH-NM com o CRB teórico para diferentes M
% ==============================================================
clear; clc;
startup;

% ------------------- Configurações ---------------------------
Mx_list = 4:1:32;           % Mx = Mz, varrendo tamanhos do array
Nreal   = 100;              % número de realizações Monte Carlo

RMSE_nm = zeros(size(Mx_list));
CRB_vec = zeros(size(Mx_list));

% -------------------- Loop por Mx = Mz ------------------------
for k = 1:numel(Mx_list)
    Mx = Mx_list(k);
    Mz = Mx;

    % Atualiza a geometria da URA
    [URA, URA_x, URA_z, x_h, x_v, z_h, z_v] = subarrays(Mx, Mz, d_x, d_z, elev, lambda, 0);
    ref = URA(1,:);  % antena de referência

    soma_nm = 0;

    % --------------- CRB para o array atual --------------------
    CRB_vec(k) = crb(L, URA, pos, lambda, P_tx, SNR_dB);

    % ---------------- Monte Carlo ------------------------------
    for r = 1:Nreal
        % Gera sinais simulados com SNR fixo
        [Yh, Yv, Y] = signals_los(pos, URA, lambda, L, alpha, SNR_dB, P_tx, Mx, Mz);

        % Estimativa inicial com PEACH
        [Un_h, Un_v, pos_est] = peach(Yh, Yv, L, ...
            x, n_hiper, x_h, z_h, x_v, z_v, ref, lambda, y, n_circ, pos);

        % Subespaço de ruído
        Cov = (Y * Y') / L;
        [V, D] = eig(Cov);
        [~, ii] = sort(diag(D), 'descend');
        Un = V(:, ii(2:end));

        % Refinamento com Nelder-Mead
        [nm_est, ~] = nelder_mead(URA, pos_est, Un, lambda, ref, ...
            deltaArea, numIterNM, 1e-6, false, x, y);

        % Erro quadrático
        soma_nm = soma_nm + norm(nm_est - pos(1:2))^2;
    end

    % RMSE final da iteração
    RMSE_nm(k) = sqrt(soma_nm / Nreal);

    
    fprintf('Mx = Mz = %d | RMSE_NM = %.4f m | CRB = %.4f m\n', Mx, RMSE_nm(k), CRB_vec(k));
end

% -------------------- Plot Final -------------------------------
figure;
semilogy(Mx_list.^2, RMSE_nm, 's-', 'LineWidth', 1.5); hold on;
semilogy(Mx_list.^2, CRB_vec, 'k--', 'LineWidth', 1.5);
xlabel('Número total de antenas M (= M_x × M_z)');
ylabel('Erro Euclidiano (m)');
grid on;
legend('PEACH-NM', 'CRB', 'Location', 'northeast');
title(sprintf('Comparação PEACH-NM vs CRB  |  M_x = M_z,  L = %d,  %d realizações, SNR = %d dB', ...
              L, Nreal, SNR_dB));
