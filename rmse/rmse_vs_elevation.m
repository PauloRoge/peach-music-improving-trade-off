clear; clc;
startup;

% ------------ VARREDURA EM ELEVAÇÃO -------------------
elev_vec = 0:1:10;  % valores de altura da URA (em metros)
RMSE = zeros(size(elev_vec));
CRBth = zeros(size(elev_vec));

% ----------------- LOOP DE ELEVAÇÕES -------------------
parfor i = 1:numel(elev_vec)
    elev = elev_vec(i);   % atualiza a altura do array

    % atualiza a posição dos elementos do array com nova altura
    [URA, ~, ~, x_h, x_v, z_h, z_v] = subarrays(Mx, Mz, ...
                                        d_x, d_z, elev, lambda, 0);

    ref = URA(1,:);  % referência permanece como a primeira antena

    CRBth(i) = crb(L, URA, pos, lambda, P_tx, SNR_dB);

    err2 = 0;

    for r = 1:MCS
        [Yh, Yv, Y] = signals_los(pos, URA, lambda, L, alpha, SNR_dB, P_tx, Mx, Mz);

        [~, ~, estPos] = peach(Yh, Yv, L, x, n_hiper, ...
            x_h, z_h, x_v, z_v, ref, lambda, y, n_circ, pos);

        %-----------------------------------------------
        % DIVISÃO DO SUBESPAÇO (Vn)
        %-----------------------------------------------
        Cov = (Y * Y') / L;
        [eigenvectors, eigenvalues] = eig(Cov); 
        estimated_sources = 1;
        [~, j] = sort(diag(eigenvalues), 'descend'); 
        eigenvectors = eigenvectors(:, j);
        Un = eigenvectors(:, estimated_sources+1:end);
        %-----------------------------------------------
        
        [nm_est, ~] = nelder_mead(URA, estPos, Un, lambda, ref, ...
            deltaArea, numIterNM, 1e-6, true, x, y);

        err2 = err2 + norm(nm_est - pos(1:2))^2;
    end

    RMSE(i) = sqrt(err2 / MCS);
    fprintf('Elev = %3d m | RMSE = %.3f m | CRB = %.3f m\n', elev, RMSE(i), CRBth(i));
end

% ---------------- PLOT RESULTADO -------------------------
figure;
semilogy(elev_vec, RMSE, '-', 'LineWidth', 1.5); hold on;
semilogy(elev_vec, CRBth, 'k--', 'LineWidth', 1.5);
xlabel('Elevação do array (m)');
ylabel('Erro (m)');
legend('RMSE Monte-Carlo','CRB teórico','Location','southwest');
grid on;
title(sprintf('PEACH-MUSIC  |  SNR = %d dB,  %d×%d URA,  L = %d,  N_{real} = %d', ...
       SNR_dB, Mx, Mz, L, MCS));