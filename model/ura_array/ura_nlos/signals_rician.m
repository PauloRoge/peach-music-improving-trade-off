% Canal Rician (LoS determinístico + difusa Rayleigh) para um único usuário.
% Mantém a mesma referência de SNR do seu código (baseline = sum(beta)).
%
% Parâmetros extras:
%   - K_dB (opcional): K-fator em dB (ex.: 8). Default = 8 dB.
%
% Chamada típica:
%   [Yh, Yv, Y] = signals_rician(UEs, URA, lambda, N, alpha, SNR_dB, P_tx, Mx, Mz, K_dB);

function [Yh, Yv, Y] = signals_rician(UEs, URA, lambda, N, ...
    alpha, SNR_dB, P_tx, Mx, Mz, K_dB)

    if nargin < 10 || isempty(K_dB)
        K_dB = 8; % LoS dominante por padrão
    end

    M           = size(URA, 1);
    user_pos    = UEs(1,:);
    SNR_linear  = 10^(SNR_dB / 10);
    K_linear    = 10^(K_dB / 10);

    % --- Geometria e perdas por elemento (mesma base do seu LoS) ---
    d_km  = sqrt(sum((URA - user_pos).^2, 2));   % [M x 1]
    d_k1  = norm(URA(1,:) - user_pos);           % ref p/ fase relativa

    % beta = (lambda / (4*pi)).^2 ./ (d_km.^alpha); % (variante potência)
    beta  = (lambda/(4*pi)).^alpha ./ (d_km.^alpha); % coerente com seu código atual

    phase = -(2*pi/lambda) * (d_k1 - d_km);
    h_los_raw = sqrt(beta) .* exp(1j * phase);   % LoS "bruto" (determinístico)
    S_baseline = sum(abs(h_los_raw).^2);         % baseline de potência (antes de K)

    % --- Alocação Rician: P_los e P_diff (difusa) ---
    P_los_target = (K_linear/(K_linear+1)) * S_baseline;
    P_dif_target = (1/(K_linear+1))       * S_baseline;

    % Escala do LoS para atingir P_los_target
    g_los = sqrt(P_los_target / max(sum(abs(h_los_raw).^2), eps));
    h_los = g_los * h_los_raw;

    % --- Componente difusa (Rayleigh por elemento, com mesma tendência de beta) ---
    % w ~ CN(0, I_M): ruído complexo branco
    w = (randn(M,1) + 1j*randn(M,1)) / sqrt(2);
    % Modelagem difusa ponderada por sqrt(beta): mantém grandezas por elemento
    h_dif_raw = sqrt(beta) .* w;
    % Normalização para atingir P_dif_target
    g_dif = sqrt(P_dif_target / max(sum(abs(h_dif_raw).^2), eps));
    h_dif = g_dif * h_dif_raw;

    % --- Canal final por elemento (Rician) ---
    H = h_los + h_dif;  % [M x 1]

    % --- Sinal transmitido e ruído (mesma SNR de referência do seu código) ---
    s = sqrt(P_tx) * (randn(1,N) + 1j*randn(1,N)) / sqrt(2);
    signal_rx = H * s;  % [M x N]

    % Referência de potência para o ruído:
    % Usa o baseline LoS puro S_baseline (mesma convenção da sua função original)
    total_rx_power_ref = P_tx * S_baseline;

    noise_power = total_rx_power_ref / SNR_linear;
    N_total     = sqrt(noise_power/2) * (randn(M,N) + 1j*randn(M,N));

    % --- Sinal final ---
    Y = signal_rx + N_total;

    % --- Extração dos subarranjos (mesma convenção) ---
    indices_h = ((1:Mx) - 1)*Mz + 1;  % primeira linha de cada coluna
    indices_v = 1:Mz;                 % primeira coluna completa

    Yh = Y(indices_h, :);
    Yv = Y(indices_v, :);
end
