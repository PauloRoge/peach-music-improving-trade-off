% Canal: LoS + L caminhos NLoS especulares (um único usuário).
% - Geometria esférica por elemento (distância euclidiana).
% - Mantém a referência de SNR do seu código original, calibrando o ruído
%   com base na potência LoS pura (baseline = sum(|h0_raw|^2)).
%
% Assinatura:
% [Yh, Yv, Y] = signals_los_plusL_specular(UEs, URA, lambda, N, ...
%     alpha, SNR_dB, P_tx, Mx, Mz, scatterers_pos, Gamma, ...
%     normalize_by_K, K_dB)
%
% Parâmetros:
% - UEs: [x y z] do usuário (usa UEs(1,:))
% - URA: [M x 3] posições dos elementos (m)
% - lambda: comprimento de onda (m)
% - N: nº snapshots
% - alpha: expoente de perda (ex.: 2)
% - SNR_dB: SNR alvo (dB)
% - P_tx: potência do símbolo s (linear)
% - Mx, Mz: dimensões do URA para extrair Yh/Yv
% - scatterers_pos: [L x 3] posições dos espalhadores (m)
% - Gamma: [L x 1] coeficientes de reflexão complexos (|Gamma|<=1)
% - normalize_by_K (opcional, default=true): habilita normalização por K
% - K_dB (opcional, default=8): K-fator LoS:(soma NLoS) em dB
%
% Saídas:
% - Y:  [M x N] sinal em todas as antenas
% - Yh: subarranjo horizontal (primeira linha de cada coluna)
% - Yv: subarranjo vertical   (primeira coluna completa)

function [Yh, Yv, Y] = signals_los_plusL_specular(UEs, URA, lambda, N, ...
    alpha, SNR_dB, P_tx, Mx, Mz, scatterers_pos, Gamma, normalize_by_K, K_dB)

    % ---------------- Defaults ----------------
    if nargin < 12 || isempty(Gamma)
        L = size(scatterers_pos,1);
        Gamma = 0.4 * exp(1j*2*pi*rand(L,1)); % reflexão moderada com fase aleatória
    end
    if nargin < 13 || isempty(normalize_by_K)
        normalize_by_K = true; % por padrão normaliza por K
    end
    if nargin < 14 || isempty(K_dB)
        K_dB = 8; % LoS dominante
    end

    % ---------------- Básico ----------------
    M          = size(URA, 1);
    user_pos   = UEs(1,:);
    SNR_linear = 10^(SNR_dB/10);
    K_linear   = 10^(K_dB/10);

    % ---------------- Caminho LoS (esférico) ----------------
    d_km_los  = sqrt(sum((URA - user_pos).^2, 2));         % [M x 1]
    d_k1_los  = norm(URA(1,:) - user_pos);                 % ref para fase
    beta_los  = (lambda/(4*pi)).^alpha ./ (d_km_los.^alpha);
    phase_los = -(2*pi/lambda) * (d_k1_los - d_km_los);
    h0_raw    = sqrt(beta_los) .* exp(1j*phase_los);       % [M x 1]

    % Baseline LoS puro para referência de SNR
    S_baseline = sum(abs(h0_raw).^2);

    % ---------------- L caminhos NLoS especulares ----------------
    L  = size(scatterers_pos, 1);
    H_spec_raw = zeros(M,1);

    for ell = 1:L
        s_pos  = scatterers_pos(ell,:);         % [1x3]
        d_us   = norm(user_pos - s_pos);        % UE -> espalhador
        d_sm   = sqrt(sum((URA - s_pos).^2, 2));% espalhador -> antena m
        d_km   = d_us + d_sm;                   % caminho total por elemento
        d_k1   = d_us + norm(s_pos - URA(1,:)); % ref para fase relativa

        beta_l = (lambda/(4*pi)).^alpha ./ (d_km.^alpha);
        phase_l= -(2*pi/lambda) * (d_k1 - d_km);

        H_spec_raw = H_spec_raw + Gamma(ell) .* sqrt(beta_l) .* exp(1j*phase_l);
    end

    % Potências "cruas" (antes de normalizações)
    P0_raw = sum(abs(h0_raw   ).^2);
    P1_raw = sum(abs(H_spec_raw).^2);

    % ---------------- Normalização por K (opcional) ----------------
    % Objetivo: definir potência LoS e soma NLoS para obedecer K,
    % mantendo a referência de SNR baseada no LoS puro (S_baseline).
    if normalize_by_K
        P0_target = (K_linear/(K_linear+1)) * S_baseline;
        P1_target = (1/(K_linear+1))       * S_baseline;

        g0 = sqrt(P0_target / max(P0_raw, eps));
        g1 = sqrt(P1_target / max(P1_raw, eps));
    else
        % Sem normalização por K: preserva |Gamma| e geometria.
        % Mantenha LoS com g0=1; ajuste g1 opcionalmente se desejar
        % que a soma (LoS+NLoS) mantenha a mesma escala típica.
        g0 = 1.0;
        g1 = 1.0;
    end

    h0     = g0 * h0_raw;        % LoS escalado
    H_spec = g1 * H_spec_raw;    % Soma dos L NLoS escalada
    H      = h0 + H_spec;        % Canal total por elemento [M x 1]

    % ---------------- Sinal transmitido e ruído ----------------
    s = sqrt(P_tx) * (randn(1,N) + 1j*randn(1,N)) / sqrt(2);
    signal_rx = H * s; % [M x N]

    % Ruído branco calibrado pela referência LoS pura (S_baseline)
    noise_power = (P_tx * S_baseline) / SNR_linear;
    N_total     = sqrt(noise_power/2) * (randn(M,N) + 1j*randn(M,N));

    % ---------------- Sinal final e subarranjos ----------------
    Y = signal_rx + N_total;

    indices_h = ((1:Mx) - 1)*Mz + 1;  % primeira linha de cada coluna
    indices_v = 1:Mz;                 % primeira coluna completa

    Yh = Y(indices_h, :);
    Yv = Y(indices_v, :);
end
