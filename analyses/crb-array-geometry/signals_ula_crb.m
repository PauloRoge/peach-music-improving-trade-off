function Y = signals_ula_crb(UEs, ULA, lambda, L, alpha, SNR_dB, P_tx, snr_mode)
% SIGNALS_ULA  Gera sinais recebidos em uma ULA (LoS, banda-estreita).
%
%   Y = signals_ula(UEs, ULA, lambda, L, alpha, SNR_dB, P_tx, snr_mode)
%
% Parâmetros:
%   UEs      : N×3 posições dos usuários (usa-se apenas o primeiro)
%   ULA      : M×3 posições dos elementos da ULA (formato [x,y,z])
%   lambda   : comprimento de onda
%   L        : nº de amostras de tempo
%   alpha    : expoente de perda de percurso (p.ex. 2 para espaço livre)
%   SNR_dB   : SNR alvo
%   P_tx     : potência média por símbolo (do sinal transmitido s)
%   snr_mode : 'per_antenna' (default) ou 'array'
%              'per_antenna' -> SNR por antena k:  P_tx*beta_k / sigma2_k = SNR
%              'array'       -> SNR agregado:      P_tx*sum(beta) / sigma2 = SNR
%
% Saída:
%   Y        : M×L sinais recebidos (canal determinístico LoS + ruído AWGN)
%
% Observações:
% - Fase de H é referenciada ao primeiro elemento da ULA (d_ref).
% - s ~ CN(0, P_tx) i.i.d.; banda-estreita (mesmo s em todas as antenas).
% - Para SNR "por antena", o ruído é gerado com variância por linha.

    if nargin < 8 || isempty(snr_mode), snr_mode = 'per_antenna'; end

    % --- Seleciona usuário e SNR ---
    user_pos   = UEs(1, :);
    SNR_linear = 10^(SNR_dB/10);

    % --- Distâncias e canal LoS ---
    % d_k: distância de cada antena ao UE; d_ref: referência (1ª antena)
    d_k   = sqrt(sum((ULA - user_pos).^2, 2));              % M×1
    d_ref = norm(ULA(1,:) - user_pos);

    % Ganho de caminho (magnitude^2) ~ (lambda/4pi)^2 / d^alpha
    c_fs  = (lambda/(4*pi))^2;
    beta  = c_fs ./ (d_k.^alpha);                           % M×1

    % Fase relativa (referência no 1º elemento)
    phase = -(2*pi/lambda) * (d_ref - d_k);                 % M×1
    H     = sqrt(beta) .* exp(1j*phase);                    % M×1

    % --- Sinal transmitido ---
    s = sqrt(P_tx) * (randn(1,L) + 1j*randn(1,L)) / sqrt(2); % 1×L, CN(0,P_tx)

    % --- Sinal recebido sem ruído ---
    signal_rx = H * s;  % (M×1)*(1×L) = M×L

    % --- Ruído AWGN conforme modo de SNR ---
    switch lower(snr_mode)
        case 'per_antenna'
            % SNR alvo por antena k:  P_tx*beta_k / sigma2_k = SNR_linear
            % => sigma2_k = P_tx*beta_k / SNR_linear
            sigma2_k = (P_tx * beta) / SNR_linear;          % M×1
            % Gerar ruído por linha com variância específica
            N = (randn(size(signal_rx)) + 1j*randn(size(signal_rx))) / sqrt(2);
            % Escala linha a linha
            N = bsxfun(@times, sqrt(sigma2_k), N);

        case 'array'
            % SNR agregado no array:  P_tx*sum(beta)/sigma2 = SNR_linear
            sigma2 = (P_tx * sum(beta)) / SNR_linear;
            N = sqrt(sigma2/2) * (randn(size(signal_rx)) + 1j*randn(size(signal_rx)));

        otherwise
            error('snr_mode deve ser ''per_antenna'' ou ''array''.');
    end

    % --- Saída final ---
    Y = signal_rx + N;
end
