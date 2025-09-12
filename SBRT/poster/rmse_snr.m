%============================================================
%  PEACH-MUSIC – Comparação de variantes com PARFOR
%============================================================
clear; clc;
start_figure2;                               % parâmetros gerais do projeto

% ---------- 1. Geometria do array ---------------------------------------
[URA, URA_x, URA_z, x_h, x_v, z_h, z_v] = ...
    subarrays(Mx, Mz, d_x, d_z, elev, lambda, 0);
ref = URA(1,:);                        % elemento de referência

% ---------- 2. Varredura de SNR -----------------------------------------
SNR_dB_vec = 0:1:20;                   % passo de 2 dB para reduzir tempo
nsnr       = numel(SNR_dB_vec);

RMSE_peach  = zeros(1,nsnr);
RMSE_nm     = zeros(1,nsnr);
RMSE_golden = zeros(1,nsnr);
CRBth       = zeros(1,nsnr);

time_pch_avg = zeros(1,nsnr);
time_nm_avg  = zeros(1,nsnr);
time_gd_avg  = zeros(1,nsnr);

time_un_avg  = zeros(1,nsnr);

% ---------- 3. Loop principal (parfor em SNR) ---------------------------
parfor k = 1:nsnr
    SNRdB = SNR_dB_vec(k);

    % acumuladores de erro e tempo para este SNR
    err2_pc = 0;
    err2_nm = 0;
    err2_hj = 0;
    err2_sb = 0;
    err2_gd = 0;
    t_pc = 0;  t_nm = 0;
    % t_hj = 0;  t_sb = 0;
    t_gd = 0;
    t_un = 0;
    crb_sum = 0;

    for r = 1:MCS
        rng(r);                                     % reprodutibilidade
        % posição aleatória do usuário
        pos = [ -50 + 100*rand , 13 + 37*rand , 0 ];

        % ---------- Sinais ------------------------------------------------
        [Yh, Yv, Y] = signals_los(pos, URA, lambda, L, ...
                              alpha, SNRdB, P_tx, Mx, Mz);

        % ---------- CRB ---------------------------------------------------
        crb_sum = crb_sum + crb(L, URA, pos, lambda, P_tx, SNRdB,2);

        % ---------- PEACH (coarse) ----------------------------------------
        tic;                                              % cronômetro
        [~, ~, est_peach] = peach( ...
            Yh, Yv, L, x, 24, x_h, z_h, x_v, z_v, ...
            ref, lambda, y, 12, pos);
        t_pc = t_pc + toc;

        % ---------- Golden-PEACH  -----------------------------
        tic;
        [~, ~, est_golden] = peach_golden( ...
            Yh, Yv, L, x, 24, x_h, z_h, x_v, z_v, ...
            ref, lambda, y, 12, pos);
        t_gd = t_gd + toc;

        % ---------- Sub-espaço de ruído para NM e Subplex ----------------
        tic;
        Cov = (Y*Y')/L;
        [V, D] = eig(Cov); [~,idx] = sort(diag(D),'descend');
        Un = V(:, idx(2:end));                          % fonte única
        t_un = t_un + toc;

        % ---------- Nelder-Mead ------------------------------------------
        tic;
        est_nm = nelder_mead(URA, est_peach, Un, lambda, ref, ...
                             deltaArea, numIterNM, tol, false, x, y);
        t_nm = t_nm + toc;

        % ---------- Acumula erros quadráticos -----------------------------
        err2_pc = err2_pc + norm(est_peach  - pos(1:2))^2;
        err2_nm = err2_nm + norm(est_nm     - pos(1:2))^2;
        err2_gd = err2_gd + norm(est_golden - pos(1:2))^2;
    end

    % ---------- Estatísticas por SNR -------------------------------------
    RMSE_peach(k)  = sqrt(err2_pc / MCS);
    RMSE_nm(k)     = sqrt(err2_nm / MCS);
    RMSE_golden(k) = sqrt(err2_gd / MCS);

    CRBth(k)       = crb_sum / MCS;

    time_un_avg(k)  = t_un / MCS;

    time_pch_avg(k) = t_pc / MCS;
    time_nm_avg(k)  = (t_nm / MCS) + time_un_avg(k) + time_pch_avg(k);
    time_gd_avg(k)  = t_gd / MCS;

    fprintf(['SNR %2d dB | RMSE→ PCH %.2f  NM %.2f' ...
             'GD %.2f | CRB %.2f  (avg de %.0f exec.)\n'], ...
            SNRdB, RMSE_peach(k), RMSE_nm(k), ...
            RMSE_golden(k), CRBth(k), MCS);
end

%% Curvas de RMSE × SNR – PEACH-MUSIC
figure('Units','centimeters','Position',[2 2 16 12]);

semilogy(SNR_dB_vec, RMSE_peach , 'x--', ...
         SNR_dB_vec, RMSE_nm    , 'o-' , ...
         SNR_dB_vec, RMSE_golden, 's--', ...
         SNR_dB_vec, CRBth      , 'k--','LineWidth',1.5);
grid on;
xlabel('SNR (dB)');
ylabel('RMSE (m)');
legend('PEACH','PEACH + Nelder-Mead','Golden PEACH','CRB','Location','best');
title(sprintf('PEACH-MUSIC – %d×%d URA  |  L=%d  |  MCS=%d', Mx,Mz,L,MCS));


% figura para artigo em duas colunas (≈8 cm cada subfigura)
fig = figure('Units','centimeters','Position',[2 2 16 7]);

tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');

% ---------- RMSE ----------
nexttile;
semilogy(SNR_dB_vec, RMSE_peach , 'x--', ...
         SNR_dB_vec, RMSE_nm    , 'o-' , ...
         SNR_dB_vec, RMSE_golden, 's--', ...
         SNR_dB_vec, CRBth      , 'k--', 'LineWidth',1.5);
grid on;
xlabel('SNR (dB)'); ylabel('RMSE (m)');
legend('PEACH','PEACH + Nelder-Mead','Golden PEACH','CRB','Location','best');
title(sprintf('PEACH-MUSIC – %dx%d URA  |  L=%d  |  MCS=%d',Mx,Mz,L,MCS));

% ---------- Tempo ----------
nexttile;
semilogy(SNR_dB_vec, time_pch_avg, 'x--', ...
         SNR_dB_vec, time_nm_avg , 'o-' , ...
         SNR_dB_vec, time_gd_avg , 's--', 'LineWidth',1.5);
grid on;
xlabel('SNR (dB)'); ylabel('Tempo médio por realização (s)');
legend('PEACH','Nelder-Mead','Golden PEACH','Location','best');
title(sprintf('Tempo médio de execução  |  %dx%d URA',Mx,Mz));

% ---------- remove espaços extra ----------
set(fig,'PaperPositionMode','auto');      % usa a janela como bounding-box
ax = findall(fig,'Type','axes');
arrayfun(@(a) set(a,'LooseInset',max(a.TightInset,0.02)), ax);  % “cola” os eixos