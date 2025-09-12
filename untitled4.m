%============================================================
%  PEACH-MUSIC – Comparação de variantes com PARFOR
%============================================================
clear; clc;
start_figure2; % parâmetros gerais do projeto

K_dB = 0;

scatterer_pos = [  % Lx3
    -35, 15, 2;
    -8,  10, 1.2;
    20, 5,  2.5
];

Gamma = [ ...
    0.5*exp(1j*pi/3);
    0.3*exp(1j*1.1*pi);
    0.4*exp(1j*0.2*pi)
];
% ---------- 1. Geometria do array ---------------------------------------
[URA, URA_x, URA_z, x_h, x_v, z_h, z_v] = ...
    subarrays(Mx, Mz, d_x, d_z, elev, lambda, 0);
ref = URA(1,:); % elemento de referência

% ---------- 2. Varredura de SNR -----------------------------------------
SNR_dB_vec = 0:1:20;                   % passo de 2 dB para reduzir tempo
nsnr       = numel(SNR_dB_vec);

%% ---------------- Pré-alocações -----------------------------------------
MSE_golden_los  = zeros(1,nsnr);
MSE_golden_nlos = zeros(1,nsnr);
CRBth           = zeros(1,nsnr);

time_gd_los_avg  = zeros(1,nsnr);
time_gd_nlos_avg = zeros(1,nsnr);

%% ---------------- Loop principal (parfor em SNR) ------------------------
parfor k = 1:nsnr
    SNRdB = SNR_dB_vec(k);

    % Acumuladores por SNR
    err2_los = 0;
    err2_nlo = 0;
    t_los    = 0;
    t_nlo    = 0;
    crb_sum  = 0;

    for r = 1:MCS
        % Semente reprodutível (única por (k,r))
        rng(10^6*k + r,'twister');

        % Posição aleatória do usuário (ajuste conforme seu cenário)
        pos = [ -50 + 100*rand , 13 + 37*rand , 0 ];

        % ---------- Sinais LoS ------------------------------------------
        % Saída esperada: Wh, Wv, Y_los (ou W) compatível com seu peach_golden
        [Wh, Wv, W_los] = signals_los(pos, URA, lambda, L, alpha, SNRdB, P_tx, Mx, Mz);

        % ---------- Sinais NLoS -----------------------------------------
        % Saída esperada: Yh, Yv, Y_nlos
        [Yh, Yv, Y_nlos] = signals_nlos_multi(pos, URA, lambda, L, alpha, SNRdB, ...
                                              P_tx, Mx, Mz, K_dB, scatterer_pos, Gamma);

        % ---------- CRB (único usuário) ---------------------------------
        crb_sum = crb_sum + crb(L, URA, pos, lambda, P_tx, SNRdB, 2);

        % ---------- PEACH-Golden (LoS) ----------------------------------
        t0 = tic;
        [~, ~, est_gd_los] = peach_golden(Yh, Yv, L, x, n_circ, x_h, z_h, x_v, z_v, ...
                                          ref, lambda, y, n_hiper, pos);
        t_los = t_los + toc(t0);

        % est_nm = nelder_mead(URA, est_pc1, Un, lambda, ref, ...
        %     deltaArea, numIterNM, tol, false, x, y);

        % ---------- PEACH-Golden (NLoS) ---------------------------------
        t0 = tic;
        [~, ~, est_gd_nlo] = peach_golden(Yh, Yv, L, x, n_circ, x_h, z_h, x_v, z_v, ...
                                          ref, lambda, y, n_hiper, pos);
        t_nlo = t_nlo + toc(t0);
        
        Cov = (Y_nlos*Y_nlos')/L;
        [V,D] = eig(Cov); [~,idx] = sort(diag(D),'descend');
        Un = V(:, idx(2:end));
        
        est_nm = nelder_mead(URA, est_gd_nlo, Un, lambda, ref, ...
            deltaArea, numIterNM, tol, false, x, y);

        % ---------- Erros quadráticos 2D (x,y) --------------------------
        err2_los = err2_los + norm(est_gd_los(1:2) - pos(1:2))^2;
        err2_nlo = err2_nlo + norm(est_nm(1:2) - pos(1:2))^2;
    end

    % Médias por SNR
    MSE_golden_los(k)  = err2_los / MCS;
    MSE_golden_nlos(k) = err2_nlo / MCS;
    CRBth(k)           = crb_sum  / MCS;

    time_gd_los_avg(k)  = t_los / MCS;
    time_gd_nlos_avg(k) = t_nlo / MCS;

    fprintf('SNR %2d dB | MSE  LoS %.3f | MSE NLoS %.3f | CRB %.3f | <t> LoS %.4fs | <t> NLoS %.4fs\n', ...
        SNRdB, MSE_golden_los(k), MSE_golden_nlos(k), CRBth(k), ...
        time_gd_los_avg(k), time_gd_nlos_avg(k));
end

%% ----------------- Curvas de MSE × SNR ----------------------------------
figure('Units','centimeters','Position',[2 2 16 12]);
semilogy(SNR_dB_vec, MSE_golden_los , 'o-', ...
         SNR_dB_vec, MSE_golden_nlos, 's-', ...
         SNR_dB_vec, CRBth          , 'k--', 'LineWidth',1.5);
grid on; xlabel('SNR (dB)'); ylabel('MSE (m^2)');
legend('Golden (LoS)','Golden (NLoS)','CRB','Location','best');
title(sprintf('PEACH-Golden – %d\\times%d URA  |  L=%d  |  MCS=%d', Mx, Mz, L, MCS));

%% ----------------- Curvas de Tempo × SNR --------------------------------
figure('Units','centimeters','Position',[2 2 16 12]);
semilogy(SNR_dB_vec, time_gd_los_avg , 'o-', ...
         SNR_dB_vec, time_gd_nlos_avg, 's-', 'LineWidth',1.5);
grid on; xlabel('SNR (dB)'); ylabel('Tempo médio por realização (s)');
legend('Golden (LoS)','Golden (NLoS)','Location','best');
title(sprintf('Tempo médio – PEACH-Golden  |  %d\\times%d URA', Mx, Mz));

% % Parâmetros de exemplo
% n_circ = [96, 48, 24];
% n_hiper   = [48, 24, 12];
% 
% % Cria a figura
% figure('Units','centimeters','Position',[2 2 16 12]);
% semilogy(SNR_dB_vec, MSE_peach , 'x-', ...
%          SNR_dB_vec, MSE_peach1, 'x--', ...
%          SNR_dB_vec, MSE_peach2, 'x--', ...
%          SNR_dB_vec, MSE_golden2,'s--',...
%          SNR_dB_vec, CRBth      , 'k--','LineWidth',1.5);
% grid on;
% xlabel('SNR (dB)');
% ylabel('MSE (m)');
% 
% % Monta as entradas da legenda com subscritos
% leg = cell(1,7);
% for i = 1:3
%     leg{i}   = sprintf('PEACH~$N_C=%d,\\;N_H=%d$', ...
%                        n_circ(i), n_hiper(i));
%     leg{i+3} = sprintf('Golden~PEACH~$N_C=%d,\\;N_H=%d$', ...
%                        n_circ(i), n_hiper(i));
% end
% leg{7} = 'CRB';
% 
% legend(leg, 'Location','best', 'Interpreter','latex');
% title(sprintf('PEACH-MUSIC – %d×%d URA  |  L=%d  |  MCS=%d', Mx, Mz, L, MCS));

% 
% % figura para artigo em duas colunas (≈8 cm cada subfigura)
% fig = figure('Units','centimeters','Position',[2 2 16 7]);
% 
% tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
% 
% % ---------- MSE ----------
% nexttile;
% semilogy(SNR_dB_vec, MSE_peach , 'x--', ...
%          SNR_dB_vec, MSE_golden, 's-', ...
%          SNR_dB_vec, MSE_golden1, 's--',...
%          SNR_dB_vec, MSE_golden2, 's--',...
%          SNR_dB_vec, CRBth      , 'k--', 'LineWidth',1.5);
% grid on;
% xlabel('SNR (dB)'); ylabel('MSE (m)');
% legend('PEACH 96-48','PEACH 48-24','PEACH 24-12','Golden PEACH 96-48', ...
%     'Golden PEACH 48-24','Golden PEACH-24-12','CRB','Location','best');
% title(sprintf('PEACH-MUSIC – %dx%d URA  |  L=%d  |  MCS=%d',Mx,Mz,L,MCS));
% 
% % ---------- Tempo ----------
% nexttile;
% semilogy(SNR_dB_vec, time_pch_avg, 'x--',...
%          SNR_dB_vec, time_pch_avg1, 'x--',...
%          SNR_dB_vec, time_pch_avg2, 'x--',...
%          SNR_dB_vec, time_gd_avg , 's--', ...
%          SNR_dB_vec, time_gd_avg1 , 's--', ...
%          SNR_dB_vec, time_gd_avg2 , 's--','LineWidth',1.5);
% grid on;
% xlabel('SNR (dB)'); ylabel('Tempo médio por realização (s)');
% legend('PEACH 96-48','PEACH 48-24','PEACH 24-12','Golden PEACH 96-48', ...
%     'Golden PEACH 48-24','Golden PEACH-24-12','Location','best');
% title(sprintf('Tempo médio de execução  |  %dx%d URA',Mx,Mz));
% 
% % ---------- remove espaços extra ----------
% set(fig,'PaperPositionMode','auto');      % usa a janela como bounding-box
% ax = findall(fig,'Type','axes');
% arrayfun(@(a) set(a,'LooseInset',max(a.TightInset,0.02)), ax);  % “cola” os eixos