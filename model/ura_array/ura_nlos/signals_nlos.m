% LoS + 1 caminho especular (NLoS) para um único usuário.
% Mantém SNR coerente com o caso LoS puro e dicionário a(x,y) baseado em LoS.
% - UEs: [x y z] do usuário (usa UEs(1,:))
% - URA: [M x 3] posições dos elementos (em metros)
% - lambda: comprimento de onda
% - N: nº de snapshots
% - alpha: expoente de perda (ex.: 2)
% - SNR_dB: SNR desejada (referência ao caso LoS puro)
% - P_tx: potência do símbolo s (linear)
% - Mx, Mz: dimensões do URA (Mx colunas, Mz linhas) para extrair Yh, Yv
% - K_dB: K-factor Rician (LoS/NLoS) em dB (ex.: 8)
% - scatterer_pos: [xs ys zs] posição do espalhador
% - Gamma: coef. de reflexão complexo (módulo <= 1), default = exp(1j*phi)
%
% Saídas: Yh, Yv, Y (M x N)
function [Yh, Yv, Y] = signals_nlos(UEs, URA, lambda, N, ...
    alpha, SNR_dB, P_tx, Mx, Mz, K_dB, scatterer_pos, Gamma)

    % ---------------- Parâmetros e entradas ----------------
    if nargin < 11 || isempty(K_dB),        K_dB = 8; end               % LoS dominante
    if nargin < 12 || isempty(scatterer_pos)
        % Exemplo: espalhador ~ parede lateral (ajuste à sua cena)
        % Coloque algo fisicamente plausível no seu cenário real.
        u_tmp = UEs(1,:);
        scatterer_pos = u_tmp + [5, 3, 1.5];  % metros
    end
    if nargin < 13 || isempty(Gamma)
        % Fase aleatória; |Gamma| moderado (0.4) por default
        Gamma = 0.4 * exp(1j*2*pi*rand);
    end

    M           = size(URA, 1);
    user_pos    = UEs(1,:);
    SNR_linear  = 10^(SNR_dB/10);
    K_linear    = 10^(K_dB/10);

    % ---------------- Caminho LoS (esférico, por elemento) ----------------
    d_km_los    = sqrt(sum((URA - user_pos).^2, 2));          % [M x 1]
    d_k1_los    = norm(URA(1,:) - user_pos);                  % ref p/ fase relativa
    beta_los    = (lambda/(4*pi)).^alpha ./ (d_km_los.^alpha);% perda por elemento
    phase_los   = -(2*pi/lambda) * (d_k1_los - d_km_los);
    h0_raw      = sqrt(beta_los) .* exp(1j*phase_los);        % [M x 1]
    P0_raw_sum  = sum(abs(h0_raw).^2);                        % soma de potências por elemento

    % ---------------- Caminho NLoS (espalhador especular) -----------------
    % Distância via espalhador: UE->S + S->elemento m
    d_us        = norm(user_pos - scatterer_pos);
    d_sm        = sqrt(sum((URA - scatterer_pos).^2, 2));     % [M x 1]
    d_km_nlos   = d_us + d_sm;                                % caminho total por elemento
    % Fase relativa: referencia no elemento 1 (consistente com o LoS)
    d_k1_nlos   = d_us + norm(scatterer_pos - URA(1,:));
    beta_nlos   = (lambda/(4*pi)).^alpha ./ (d_km_nlos.^alpha);
    phase_nlos  = -(2*pi/lambda) * (d_k1_nlos - d_km_nlos);
    h1_raw      = Gamma .* sqrt(beta_nlos) .* exp(1j*phase_nlos);  % inclui reflexão
    P1_raw_sum  = sum(abs(h1_raw).^2);

    % ---------------- Particionamento Rician (K-factor) -------------------
    % Objetivo: manter SNR coerente com o seu cálculo original (caso LoS puro).
    % Usamos como "baseline" a potência LoS pura: S_baseline = sum(|h0_raw|^2).
    % Alocamos P0_target = K/(K+1)*S_baseline e P1_target = 1/(K+1)*S_baseline.
    S_baseline  = P0_raw_sum;
    P0_target   = (K_linear/(K_linear+1)) * S_baseline;
    P1_target   = (1/(K_linear+1))       * S_baseline;

    % Ganhos de normalização por caminho para atingir o K e manter baseline
    g0 = sqrt(P0_target / max(P0_raw_sum, eps));
    g1 = sqrt(P1_target / max(P1_raw_sum, eps));

    h0 = g0 * h0_raw;
    h1 = g1 * h1_raw;

    % Canal final por elemento (LoS + NLoS)
    H  = h0 + h1;                           % [M x 1]

    % ---------------- Sinal transmitido e ruído ---------------------------
    s = sqrt(P_tx) * (randn(1,N) + 1j*randn(1,N)) / sqrt(2);  % [1 x N]
    signal_rx     = H * s;                                    % [M x N]

    % Potência total-alvo para cálculo do ruído:
    % Mantemos o mesmo baseline de potência do caso LoS puro:
    total_rx_power_ref = P_tx * S_baseline;

    noise_power = total_rx_power_ref / SNR_linear;
    N_total     = sqrt(noise_power/2) * (randn(M,N) + 1j*randn(M,N));

    % ---------------- Sinal final e subarranjos ---------------------------
    Y = signal_rx + N_total;

    % Extração dos sub-arranjos (mesma convenção da sua função original)
    indices_h = ((1:Mx) - 1)*Mz + 1;   % primeira linha de cada coluna
    indices_v = 1:Mz;                  % primeira coluna completa

    Yh = Y(indices_h, :);
    Yv = Y(indices_v, :);
end
