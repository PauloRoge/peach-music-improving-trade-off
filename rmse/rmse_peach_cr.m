% ==============================================================
%  -– compara PEACH-MUSIC com o CRB teórico corretamente
% ==============================================================
clear;  clc;
startup;

% -------------- 2. PREPARA ARRAY (como seu subarrays.m) --------
[URA, URA_x, URA_z, x_h, x_v, z_h, z_v] = subarrays(Mx, Mz, d_x, d_z, elev, lambda, 0);

ref = URA(1,:);   % referência = primeiro elemento

% -------------- 3. LAÇO DE SNRs --------------------------------
SNR_dB_vec = 0:1:20;            % valores de SNR a testar
RMSE  = zeros(size(SNR_dB_vec));
CRBth = zeros(size(SNR_dB_vec));

% --- Monte-Carlo ---
for k = 1:numel(SNR_dB_vec)
    SNRdB = SNR_dB_vec(k); 
    err2  = 0;                       % acumula erro²

    % ------ CRB para este SNR ------
    CRBth(k) = crb(L, URA, pos, lambda, P_tx, SNRdB,alpha);

    % ------ realizações ------------
    for r = 1:MCS
        % gera sinal com SNR global (signals_snr.m)
        [Yh,Yv,Y] = signals(pos,URA,lambda,L,alpha,...
                                SNRdB,P_tx,Mx,Mz);

        % -------- estimador PEACH --------
        [~,~,estPeach] = peach_analitico(Yh,Yv,L, ...
            x, n_hiper, ...
            x_h, z_h, x_v, z_v, ref, ...
            lambda, y, n_circ, pos);

        %-----------------------------------------------
        % DIVISÃO DO SUBESPAÇO (Vn)
        %-----------------------------------------------
        Cov = (Y * Y') / L;
        [eigenvectors, eigenvalues] = eig(Cov); 
        estimated_sources = 1;
        [~, i] = sort(diag(eigenvalues), 'descend'); 
        eigenvectors = eigenvectors(:, i);
        Un = eigenvectors(:, estimated_sources+1:end);
        %-----------------------------------------------

        % [nm_est, simplex_history] = nelder_mead(URA, estPeach, Un, lambda, ref, ...
        % deltaArea, numIterNM, 1e-6, true, x, y);

        err2 = err2 + norm(estPeach - pos(1:2))^2;
    end

    RMSE(k) = sqrt(err2/MCS);
    fprintf('SNR = %+3d dB | RMSE = %.3f m | CRB = %.3f m\n',...
             SNRdB, RMSE(k), CRBth(k));
end

% -------------- 4. PLOT ----------------------------------------
figure;  semilogy(SNR_dB_vec,RMSE,'o-','LineWidth',1.5); hold on;
semilogy(SNR_dB_vec,CRBth,'s--','LineWidth',1.5);
grid on; xlabel('SNR (dB)'); ylabel('Erro (m)');
legend('RMSE Monte-Carlo','CRB teórico','Location','southwest');
title(sprintf('PEACH-MUSIC  |  %d×%d URA,  L = %d,  N_{real} = %d',...
               Mx,Mz,L,MCS));