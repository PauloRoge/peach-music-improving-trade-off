% LoS + L caminhos especulares (NLoS) para um único usuário.
% - Aceita múltiplos espalhadores: scatterer_pos é [L x 3], Gamma é [L x 1] (ou escalar).
% - Mantém SNR coerente com o caso LoS puro (baseline = sum(|h0_raw|^2)).
%
% Assinatura:
% [Yh, Yv, Y] = signals_nlos_multi(UEs, URA, lambda, N, ...
%   alpha, SNR_dB, P_tx, Mx, Mz, K_dB, scatterer_pos, Gamma)
%
% Parâmetros:
% - K_dB: K-factor Rician (LoS/NLoS) em dB (ex.: 8)
% - scatterer_pos: [L x 3] posições dos espalhadores (cada linha: [xs ys zs])
% - Gamma: [L x 1] coeficientes de reflexão complexos (ou escalar)
%
% Saídas: Yh, Yv, Y (M x N)

function [Yh, Yv, Y] = signals_nlos_multi(UEs, URA, lambda, N, ...
    alpha, SNR_dB, P_tx, Mx, Mz, K_dB, scatterer_pos, Gamma)

    % ---------------- Parâmetros e defaults ----------------
    if nargin < 10 || isempty(K_dB),        K_dB = 8; end
    M           = size(URA, 1);
    user_pos    = UEs(1,:);
    SNR_linear  = 10^(SNR_dB/10);
    K_linear    = 10^(K_dB/10);

    % Se não vier espalhador, cria 1 default
    if nargin < 11 || isempty(scatterer_pos)
        u_tmp = user_pos;
        scatterer_pos = u_tmp + [5, 3, 1.5]; % [1 x 3]
    end
    % Garante matriz Lx3
    if size(scatterer_pos,2) ~= 3
        error('scatterer_pos deve ser Lx3 (cada linha: [xs ys zs]).');
    end
    L = size(scatterer_pos,1);

    % Gamma pode ser escalar ou Lx1; se vazio, define aleatório moderado
    if nargin < 12 || isempty(Gamma)
        Gamma = 0.4 * exp(1j*2*pi*rand(L,1));
    elseif isscalar(Gamma)
        Gamma = repmat(Gamma, L, 1);
    else
        if numel(Gamma) ~= L
            error('Gamma deve ser escalar ou vetor Lx1 para L espalhadores.');
        end
        Gamma = Gamma(:);
    end

    % ---------------- Caminho LoS (esférico, por elemento) ----------------
    d_km_los   = sqrt(sum((URA - user_pos).^2, 2));   % [M x 1]
    d_k1_los   = norm(URA(1,:) - user_pos);           % ref p/ fase
    beta_los = (lambda / (4*pi)).^2 ./ (d_km_los.^alpha);
    phase_los  = -(2*pi/lambda) * (d_k1_los - d_km_los);
    h0_raw     = sqrt(beta_los) .* exp(1j*phase_los); % [M x 1]
    S_baseline = sum(abs(h0_raw).^2);                 % baseline LoS puro

    % ---------------- Caminhos NLoS (L especulares) -----------------------
    hN_raw = zeros(M, L);
    Pl_raw = zeros(L,1);
    for l = 1:L
        s_l      = scatterer_pos(l,:);
        d_us     = norm(user_pos - s_l);                 % escalar
        d_sm     = sqrt(sum((URA - s_l).^2, 2));         % [M x 1]
        d_ml     = d_us + d_sm;                          % [M x 1]
        d_ref_l  = d_us + norm(s_l - URA(1,:));          % ref no elemento 1
        beta_l   = (lambda/(4*pi)).^alpha ./ (d_ml.^alpha);
        phase_l  = -(2*pi/lambda) * (d_ref_l - d_ml);
        h_l_raw  = Gamma(l) .* sqrt(beta_l) .* exp(1j*phase_l); % [M x 1]
        hN_raw(:,l) = h_l_raw;
        Pl_raw(l)   = sum(abs(h_l_raw).^2);              % potência "bruta" do caminho l
    end

    % ---------------- Particionamento Rician (K-factor) -------------------
    P0_target  = (K_linear/(K_linear+1)) * S_baseline;
    PN_target  = (1/(K_linear+1))       * S_baseline;

    % Normaliza LoS
    g0   = sqrt(P0_target / max(sum(abs(h0_raw).^2), eps));
    h0   = g0 * h0_raw;

    % Distribui PN_target entre os L caminhos proporcionalmente a Pl_raw
    if all(Pl_raw <= 0)
        % fallback: distribuição uniforme
        wl = ones(L,1) / L;
    else
        wl = Pl_raw / sum(Pl_raw);  % pesos somam 1
    end

    % Normaliza cada caminho NLoS para atingir PN_target*wl(l)
    hN = zeros(M,1);
    for l = 1:L
        target_l = PN_target * wl(l);
        gl = sqrt(target_l / max(Pl_raw(l), eps));
        hN = hN + gl * hN_raw(:,l);
    end

    % Canal final por elemento
    H = h0 + hN;   % [M x 1]

    % ---------------- Sinal transmitido e ruído ---------------------------
    s = sqrt(P_tx) * (randn(1,N) + 1j*randn(1,N)) / sqrt(2);  % [1 x N]
    signal_rx = H * s;                                        % [M x N]

    % Ruído calibrado com baseline LoS puro (mantém SNR compatível)
    total_rx_power_ref = P_tx * S_baseline;
    noise_power = total_rx_power_ref / SNR_linear;
    N_total     = sqrt(noise_power/2) * (randn(M,N) + 1j*randn(M,N));

    % ---------------- Sinal final e subarranjos ---------------------------
    Y = signal_rx + N_total;

    indices_h = ((1:Mx) - 1)*Mz + 1;   % primeira linha de cada coluna
    indices_v = 1:Mz;                  % primeira coluna completa

    Yh = Y(indices_h, :);
    Yv = Y(indices_v, :);
end
