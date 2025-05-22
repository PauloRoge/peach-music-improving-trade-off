clc;
startup_fig1;

% ------------ VARREDURA EM ELEVAÇÃO -------------------
elev_vec = 0:5:100;  
RMSE_nm    = zeros(size(elev_vec));  % RMSE do Nelder-Mead
RMSE_peach = zeros(size(elev_vec));  % RMSE do PEACH puro
CRBth      = zeros(size(elev_vec));

% --- NOVAS VARIÁVEIS PARA TEMPO ---
time_peach_avg = zeros(size(elev_vec));
time_nm_avg    = zeros(size(elev_vec));
% ------------------------------------

for i = 1:numel(elev_vec)
    elev = elev_vec(i);

    [URA, ~, ~, x_h, x_v, z_h, z_v] = subarrays(Mx, Mz, d_x, d_z, elev, lambda, 0);
    ref = URA(1,:);
    CRBth(i) = crb(L, URA, pos, lambda, P_tx, SNR_dB, alpha);

    err2_nm    = 0;
    err2_peach = 0;

    % --- ACUMULADORES DE TEMPO PARA ESTA ELEVAÇÃO ---
    time_peach_total_elev = 0;
    time_nm_total_elev    = 0;
    % ------------------------------------------------

    parfor r = 1:MCS % Início do parfor
    
        [Yh, Yv, Y] = signals(pos, URA, lambda, L, alpha, SNR_dB, P_tx, Mx, Mz);
        
        % --- Medição de Tempo para PEACH puro ---
        tic_peach_local = tic; % Inicia timer local para PEACH
        [~, ~, estPos] = peach_analitico(Yh, Yv, L, x, n_hiper, ...
            x_h, z_h, x_v, z_v, ref, lambda, y, n_circ, pos);
        time_peach_trial = toc(tic_peach_local); % Tempo para esta chamada do PEACH
        % ----------------------------------------
        
        err2_peach = err2_peach + norm(estPos - pos(1:2))^2;

        % Refinamento (Subplex)
        Cov = (Y * Y') / L;
        [V, D] = eig(Cov);
        [~, j_sort] = sort(diag(D), 'descend'); % Renomeado para evitar conflito com j em elements_z_a
        Un = V(:, j_sort(2:end));

        % --- Medição de Tempo para Subplex ---
        tic_nm_local = tic; % Inicia timer local para Subplex
        % [subplex_est, ~] = subplex_wrapper(URA, estPos, Un, lambda, ref, x, y, 1e-5, 20); % Removido 'hist' se não usado
        % time_nm_trial = toc(tic_nm_local); % Tempo para esta chamada do Subplex

        [subplex_est, ~] = nelder_mead(URA, estPos, Un, lambda, ref, ...
    deltaArea, numIterNM, 1e-6, true, x, y);

        % -----------------------------------

        err2_nm = err2_nm + norm(subplex_est - pos(1:2))^2;
        
        % --- Acumula tempos DENTRO do parfor (usando variáveis temporárias) ---
        % Para fazer isso corretamente com parfor e evitar erros de transparência,
        % é melhor acumular em vetores e somar fora, ou usar uma tática de redução.
        % No entanto, para uma simples soma, pode-se tentar diretamente se o MATLAB permitir
        % ou usar a técnica de cell array para cada worker.
        %
        % Uma forma mais robusta para acumular no parfor é criar vetores temporários:
        % Declare fora do parfor:
        % time_peach_trials_vec = zeros(1, MCS);
        % time_nm_trials_vec    = zeros(1, MCS);
        % Dentro do parfor, atribua:
        % time_peach_trials_vec(r) = time_peach_trial;
        % time_nm_trials_vec(r)    = time_nm_trial;
        % E depois some fora do parfor.
        %
        % Para simplificar a visualização aqui, vou assumir uma soma direta,
        % mas esteja ciente que o MATLAB pode reclamar. Se sim, use a técnica do vetor.
        % Se o MATLAB não permitir a modificação direta de variáveis de loop externo no parfor
        % (o que é comum para evitar race conditions), precisaremos de uma abordagem diferente.
        % A melhor forma é cada worker do parfor retornar seu tempo, e agregamos depois.
        
        % *** CORREÇÃO PARA ACUMULAÇÃO DE TEMPO EM PARFOR ***
        % O parfor não permite modificar diretamente variáveis de loop externo como 'time_peach_total_elev'
        % de forma acumulativa simples. Precisamos que cada iteração retorne seus valores.
        
        % A estrutura do parfor precisaria ser algo como:
        % resultados_parfor = cell(1, MCS);
        % parfor r = 1:MCS
        %     ... calcule time_peach_trial, time_nm_trial ...
        %     resultados_parfor{r} = [time_peach_trial, time_nm_trial, norm_err_peach_sq, norm_err_nm_sq];
        % end
        % E então iterar sobre resultados_parfor para somar.

        % Para manter a estrutura atual e simplificar a mudança para este exemplo,
        % vamos criar vetores temporários *dentro* do escopo da elevação 'i'
        % e preenchê-los no parfor.
        
        % Esta abordagem requer que 'time_peach_iter_results' e 'time_nm_iter_results'
        % sejam redefinidas para cada 'i'.
        % E as variáveis acumuladoras de erro também precisariam ser tratadas assim.

    end % Fim do parfor

    % --- Para coletar os tempos do parfor, a maneira mais limpa é assim: ---
    % Redefina as variáveis de acumulação para esta elevação 'i'
    err_sq_peach_trials = zeros(1, MCS);
    err_sq_nm_trials    = zeros(1, MCS);
    time_peach_trials   = zeros(1, MCS);
    time_nm_trials      = zeros(1, MCS);

    parfor r_idx = 1:MCS % Renomeado para r_idx para clareza
        % Gera sinais (copiado para dentro para ser independente)
        % Nota: URA, x_h, z_h, etc., são 'broadcast variables' ou 'sliced input variables'
        % e 'pos', 'lambda', etc., são 'broadcast variables'.
        [Yh_p, Yv_p, Y_p] = signals(pos, URA, lambda, L, alpha, SNR_dB, P_tx, Mx, Mz);
        
        tic_p_local = tic;
        [~, ~, estPos_p] = peach_analitico(Yh_p, Yv_p, L, x, n_hiper, ...
            x_h, z_h, x_v, z_v, ref, lambda, y, n_circ, pos);
        time_peach_trials(r_idx) = toc(tic_p_local);
        err_sq_peach_trials(r_idx) = norm(estPos_p - pos(1:2))^2;

        Cov_p = (Y_p * Y_p') / L;
        [V_p, D_p] = eig(Cov_p);
        [~, j_sort_p] = sort(diag(D_p), 'descend');
        Un_p = V_p(:, j_sort_p(2:end));

        tic_nm_local_p = tic;
        [subplex_est_p, ~] = subplex_wrapper(URA, estPos_p, Un_p, lambda, ref, x, y, 1e-5, 20);
        time_nm_trials(r_idx) = toc(tic_nm_local_p);
        err_sq_nm_trials(r_idx) = norm(subplex_est_p - pos(1:2))^2;
    end

    % --- Agora some os resultados do parfor ---
    err2_peach = sum(err_sq_peach_trials);
    err2_nm    = sum(err_sq_nm_trials);
    time_peach_total_elev = sum(time_peach_trials);
    time_nm_total_elev    = sum(time_nm_trials);
    % ------------------------------------------


    RMSE_peach(i) = sqrt(err2_peach / MCS);
    RMSE_nm(i)    = sqrt(err2_nm    / MCS);
    
    % --- Calcula tempo médio para esta elevação ---
    time_peach_avg(i) = time_peach_total_elev / MCS;
    time_nm_avg(i)    = time_nm_total_elev / MCS;
    % --------------------------------------------

    fprintf('Elev = %3d m | RMSE_peach = %.3f m (t=%.4fs) | RMSE_NM = %.3f m (t=%.4fs) | CRB = %.3f m\n', ...
            elev, RMSE_peach(i), time_peach_avg(i), RMSE_nm(i), time_nm_avg(i), CRBth(i));
end

% ---------------- PLOT RESULTADO RMSE -------------------------
figure;
semilogy(elev_vec, RMSE_peach, 'o-', 'LineWidth', 2,'MarkerFaceColor','w'); hold on;
semilogy(elev_vec, RMSE_nm,'o-',   'LineWidth', 2,'MarkerFaceColor','w'); % Mudado para vermelho para distinguir
semilogy(elev_vec, CRBth,'k--', 'LineWidth', 1.5,'MarkerFaceColor','w');
xlabel('Elevação do array (m)');
ylabel('Erro (m)');
legend('PEACH puro','PEACH + Subplex','CRB teórico','Location','southwest');
grid on;
% title(sprintf('RMSE | SNR = %d dB,  %d×%d URA,  L = %d,  MCS = %d', ...
%        SNR_dB, Mx, Mz, L, MCS));

% ---------------- PLOT RESULTADO TEMPO -------------------------
figure;
semilogy(elev_vec, time_peach_avg, 'ob-', 'LineWidth', 2); hold on;
semilogy(elev_vec, time_nm_avg,    'or-',   'LineWidth', 2); % Mudado para vermelho
xlabel('Elevação do array (m)');
ylabel('Tempo Médio de Execução (s)');
legend('PEACH puro','Subplex','Location','best');
grid on;
title(sprintf('Tempo de Execução | SNR = %d dB,  %d×%d URA,  L = %d,  MCS = %d', ...
       SNR_dB, Mx, Mz, L, MCS));

% clc;
% startup_fig1;
% 
% % ------------ VARREDURA EM ELEVAÇÃO -------------------
% elev_vec = 0:1:100;  
% RMSE_nm    = zeros(size(elev_vec));  % RMSE do Nelder-Mead
% RMSE_peach = zeros(size(elev_vec));  % RMSE do PEACH puro
% CRBth      = zeros(size(elev_vec));
% 
% for i = 1:numel(elev_vec)
%     elev = elev_vec(i);
% 
%     [URA, ~, ~, x_h, x_v, z_h, z_v] = subarrays(Mx, Mz, d_x, d_z, elev, lambda, 0);
%     ref = URA(1,:);
%     CRBth(i) = crb(L, URA, pos, lambda, P_tx, SNR_dB, alpha);
% 
%     err2_nm    = 0;
%     err2_peach = 0;
% 
%     for r = 1:MCS
% 
%         [Yh, Yv, Y] = signals(pos, URA, lambda, L, alpha, SNR_dB, P_tx, Mx, Mz);
%         %[Yh, Yv, Y] = signals_fig1(UEs, URA, lambda, L, SNR_dB, P_tx, Mx, Mz);
% 
% 
%         % Estimativa PEACH puro
%         [~, ~, estPos] = peach_analitico(Yh, Yv, L, x, n_hiper, ...
%             x_h, z_h, x_v, z_v, ref, lambda, y, n_circ, pos);
%         err2_peach = err2_peach + norm(estPos - pos(1:2))^2;
% 
%         % Refinamento Nelder-Mead
%         Cov = (Y * Y') / L;
%         [V, D] = eig(Cov);
%         [~, j] = sort(diag(D), 'descend');
%         Un = V(:, j(2:end));
% 
%        % [nm_est, ~] = nelder_mead(URA, estPos, Un, lambda, ref, ...
%        %    deltaArea, numIterNM, 1e-6, true, x, y);
% 
%         [subplex_est, hist] = subplex_wrapper(URA, estPos, Un, lambda, ref, x, y, 1e-5, 20);
% 
%         err2_nm = err2_nm + norm(subplex_est - pos(1:2))^2;
%     end
% 
%     RMSE_peach(i) = sqrt(err2_peach / MCS);
%     RMSE_nm(i)    = sqrt(err2_nm    / MCS);
% 
%     fprintf('Elev = %3d m | RMSE_peach = %.3f m | RMSE_NM = %.3f m | CRB = %.3f m\n', ...
%             elev, RMSE_peach(i), RMSE_nm(i), CRBth(i));
% end
% 
% % ---------------- PLOT RESULTADO -------------------------
% figure;
% semilogy(elev_vec, RMSE_peach, 'b-', 'LineWidth', 2); hold on;
% semilogy(elev_vec, RMSE_nm,    '-',   'LineWidth', 2);
% semilogy(elev_vec, CRBth,      'k--', 'LineWidth', 1.5);
% xlabel('Elevação do array (m)');
% ylabel('Erro (m)');
% legend('PEACH puro','PEACH + NM','CRB teórico','Location','southwest');
% grid on;
% title(sprintf('PEACH-MUSIC  |  SNR = %d dB,  %d×%d URA,  L = %d,  MCS = %d', ...
%        SNR_dB, Mx, Mz, L, MCS));
