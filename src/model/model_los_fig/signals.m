% essa function utiliza o mesmo metodo para calcular a snr do script do
% Bruno + a correção da antena referencia.
function [Yh, Yv, Y] = signals(UEs, URA, lambda, L, alpha, SNR_dB, P_tx, Mx, Mz)

    M = size(URA, 1);
    user_pos = UEs(1,:);                         
    SNR_linear = 10^(SNR_dB / 10);

    % --- Canal vetorizado (H) ---
    d_km  = sqrt(sum((URA - user_pos).^2, 2));
    d_k1  = norm(URA(1,:) - user_pos);

    %beta  = (lambda^alpha) ./ ((4*pi*d_km).^alpha);
    beta = (lambda / (4*pi)).^2 ./ (d_km.^alpha); % potência
    phase = -(2*pi/lambda) * (d_k1 - d_km);
    H     = sqrt(beta) .* exp(1j * phase);       

    % --- Sinal transmitido ---
    s = sqrt(P_tx) * (randn(1,L) + 1j*randn(1,L)) / sqrt(2);

    % --- Sinal recebido ideal ---
    signal_rx = H * s;                           

    % --- Cálculo de potência ---
    total_rx_power = P_tx * sum(beta);           

    % --- Ruído branco aditivo ---
    noise_power = total_rx_power / SNR_linear;
    N_total    = sqrt(noise_power/2) * (randn(M,L) + 1j*randn(M,L));

    % --- Sinal final com ruído ---
    Y = signal_rx + N_total;

    % --- Extração dos sub-arranjos ---
    indices_h = ((1:Mx) - 1)*Mz + 1;    
    indices_v = 1:Mz;                  

    Yh = Y(indices_h, :);              
    Yv = Y(indices_v, :);              

end
