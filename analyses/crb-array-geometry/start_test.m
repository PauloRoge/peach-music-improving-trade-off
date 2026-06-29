lambda = 0.03;
ULA    = ula(16, 0.5*lambda, 2.0, lambda, false);   % M×3
UEs    = [30, 0, 2.0];                              % um usuário
L      = 200; alpha = 2; SNR_dB = 10; P_tx = 1;

% SNR por antena (recomendado, consistente com seu cabeçalho antigo)
Y = signals_ula_crb(UEs, ULA, lambda, L, alpha, SNR_dB, P_tx, 'per_antenna');
