% ==========================================================
% signals.m
% Gera os sinais recebidos num cenário LoS banda-estreita.
%
% Uso:
%   [Yh, Yv, Y] = signals(UEs, URA, ULA_h, ULA_v, ...
%                         lambda, L, alpha, SNR_dB, ...
%                         P_tx, Mx, Mz, model)
%
% - UEs      : N×3 posições dos usuários (usa-se apenas o 1º)
% - URA      : M×3 posições do array completo (32×32, por ex.)
% - ULA_h    : Mx×3 posições do ULA horizontal (saída do subarrays_ula)
% - ULA_v    : Mz×3 posições do ULA vertical   (saída do subarrays_ula)
% - lambda   : comprimento de onda
% - L        : nº de amostras de tempo
% - alpha    : expoente de perda de percurso
% - SNR_dB   : SNR alvo (por antena, LoS)
% - P_tx     : potência média transmitida (por símbolo)
% - Mx,Mz    : nº de elementos nos ULAs h e v
% - model    : 'URA'  (padrão) ou 'Lshape'
%
% Saída:
%   Yh  – sinais nas antenas do ULA horizontal
%   Yv  – sinais nas antenas do ULA vertical
%   Y   – matriz completa (URAs ou ULAs empilhadas)
% ==========================================================
function [Yh, Yv, Y] = signals_ula(UEs, URA, ULA_h, ULA_v, ...
                                lambda, L, alpha, SNR_dB, ...
                                P_tx, Mx, Mz, model)

    if nargin < 13 || isempty(model)
        model = 'URA';
    end

    user_pos   = UEs(1, :);
    SNR_linear = 10^(SNR_dB/10);

    switch lower(model)

        % --------------------------------------------------
        % MODELO PARA URA COMPLETO (e.g. 32×32 = 1024)
        % --------------------------------------------------
        case 'ura'

            M = size(URA, 1);

            d_km  = sqrt(sum((URA - user_pos).^2, 2));
            d_ref = norm(URA(1,:) - user_pos);

            beta  = (lambda/(4*pi)).^2 ./ (d_km.^alpha);
            phase = -(2*pi/lambda) * (d_ref - d_km);
            H     = sqrt(beta) .* exp(1j*phase);

            s          = sqrt(P_tx)*(randn(1,L) + 1j*randn(1,L))/sqrt(2);
            signal_rx  = H * s;

            total_rx_power = P_tx * sum(beta);
            noise_power    = total_rx_power / SNR_linear;
            N_total        = sqrt(noise_power/2)*(randn(M,L)+1j*randn(M,L));

            Y = signal_rx + N_total;

            % índices dos sub-arranjos dentro do URA
            idx_h = ((1:Mx) - 1)*Mz + 1;   % 1ª coluna de cada linha
            idx_v = 1:Mz;                  % 1ª linha do URA

            Yh = Y(idx_h, :);
            Yv = Y(idx_v, :);

        % --------------------------------------------------
        % MODELO PARA L-SHAPE (ULA_h + ULA_v),  Mx+Mz antenas
        % --------------------------------------------------
        case 'lshape'

            % ---------- ULA horizontal ----------
            d_h   = sqrt(sum((ULA_h - user_pos).^2, 2));
            d_ref = norm(ULA_h(1,:) - user_pos);             % canto comum
            beta_h  = (lambda/(4*pi)).^2 ./ (d_h.^alpha);
            phase_h = -(2*pi/lambda)*(d_ref - d_h);
            H_h     = sqrt(beta_h).*exp(1j*phase_h);

            % ---------- ULA vertical -------------
            d_v   = sqrt(sum((ULA_v - user_pos).^2, 2));
            beta_v  = (lambda/(4*pi)).^2 ./ (d_v.^alpha);
            phase_v = -(2*pi/lambda)*(d_ref - d_v);
            H_v     = sqrt(beta_v).*exp(1j*phase_v);

            % ---------- símbolo transmitido -----
            s = sqrt(P_tx)*(randn(1,L) + 1j*randn(1,L))/sqrt(2);

            signal_rx_h = H_h * s;
            signal_rx_v = H_v * s;

            total_rx_power = P_tx*(sum(beta_h) + sum(beta_v));
            noise_power    = total_rx_power / SNR_linear;

            N_h = sqrt(noise_power/2)*(randn(Mx,L)+1j*randn(Mx,L));
            N_v = sqrt(noise_power/2)*(randn(Mz,L)+1j*randn(Mz,L));

            Yh = signal_rx_h + N_h;
            Yv = signal_rx_v + N_v;

            % empilha para ter a “visão” completa do array
            Y  = [Yh; Yv];

        otherwise
            error('Modelo desconhecido: use ''URA'' ou ''Lshape''.');
    end
end
