function [crb_euc, crb_x, crb_y, FIM] = crb_nlos(L, URA, user_pos, lambda, P_tx, SNR_dB, alpha, K_dB, scatterer_pos, Gamma)
% CRB para posição 2D (x,y) no modelo LoS + L especulares (NLoS), 
% consistente com signals_nlos_multi.
%
% Assinatura:
%   [crb_euc, crb_x, crb_y, FIM] = crb_nlos(L, URA, user_pos, lambda, P_tx, SNR_dB, alpha, K_dB, scatterer_pos, Gamma)
%
% Parâmetros:
%   L               : nº de snapshots
%   URA [M x 3]     : coordenadas dos elementos do array
%   user_pos [1x3]  : posição do usuário [x y z]
%   lambda          : comprimento de onda
%   P_tx            : potência por snapshot (consistente com seu 's')
%   SNR_dB          : SNR em dB (definida como no signals_nlos_multi)
%   alpha           : expoente de perda de percurso
%   K_dB            : fator Rician (dB)
%   scatterer_pos   : [Lscat x 3] posições dos espalhadores
%   Gamma           : [Lscat x 1] coeficientes complexos dos espalhadores (ou escalar)
%
% Saídas:
%   crb_euc         : CRB_x + CRB_y (m^2)
%   crb_x, crb_y    : termos marginais
%   FIM             : matriz de informação de Fisher 2x2 em (x,y)
%
% Notas:
% - Este CRB inclui as derivadas de g_l (normalização NLoS), evitando viés
%   por ignorar o reescalonamento de potência entre os caminhos.
% - Segue a calibração do ruído do seu gerador: sigma^2 = P_tx*S_baseline/SNR_lin.
% - Para refletir exatamente seu gerador:
%   * LoS usa beta_los = (lambda/(4*pi))^2 / d^alpha
%   * NLoS usa beta_l  = (lambda/(4*pi))^alpha / d^alpha  (como no arquivo)
%   Se preferir consistência física (^(2) em ambos), troque 'Calpha' para 'C2'.
%
% Referência de implementação do gerador: signals_nlos_multi  (mantida a coerência de SNR/partição K).

    % -------------------- Preparos e constantes --------------------------
    M        = size(URA,1);
    x        = user_pos(1);  y = user_pos(2);  z = user_pos(3);
    x_ref    = URA(1,1);     y_ref = URA(1,2); z_ref = URA(1,3);

    SNR_lin  = 10.^(SNR_dB/10);
    K_lin    = 10.^(K_dB/10);

    C2       = (lambda/(4*pi))^2;     % constante LoS (potência)
    Calpha   = (lambda/(4*pi))^alpha; % constante NLoS (potência) - espelha seu gerador

    % Vetores auxiliares
    dxm = x - URA(:,1);
    dym = y - URA(:,2);
    dzm = z - URA(:,3);
    dm  = sqrt(dxm.^2 + dym.^2 + dzm.^2) + eps; % distâncias usuário->antena m
    dref = sqrt((x - x_ref)^2 + (y - y_ref)^2 + (z - z_ref)^2) + eps;

    % -------------------- LoS "raw" (antes da partição K) ----------------
    beta0_raw = C2 ./ (dm.^alpha);        % |h0_raw|^2 por elemento
    a0_raw    = sqrt(beta0_raw);          % amplitude
    phi0      = (2*pi/lambda) * (dm - dref);
    h0_raw    = a0_raw .* exp(1j*phi0);   % [M x 1]

    % Potência de referência LoS (baseline) e ruído (como no gerador)
    S_baseline    = sum(abs(h0_raw).^2);                 % = sum(beta0_raw)
    noise_power   = (P_tx * S_baseline) ./ SNR_lin;      % sigma^2 por amostra
    Es            = L * P_tx;                            % energia total do sinal conhecido
    pre_factor    = (2 * Es) ./ noise_power;             % fator da FIM (Gaussiano complexo)

    % -------------------- Caminhos NLoS "raw" ----------------------------
    Lscat = size(scatterer_pos,1);
    if isscalar(Gamma), Gamma = repmat(Gamma, Lscat, 1); else, Gamma = Gamma(:); end

    % Distâncias usuário->espalhador (dependem do usuário)
    dus  = sqrt( (x - scatterer_pos(:,1)).^2 + ...
                 (y - scatterer_pos(:,2)).^2 + ...
                 (z - scatterer_pos(:,3)).^2 ) + eps;          % [Lscat x 1]
    % Derivadas de d_us
    ddus_dx = (x - scatterer_pos(:,1)) ./ dus;                  % [Lscat x 1]
    ddus_dy = (y - scatterer_pos(:,2)) ./ dus;

    % Distâncias espalhador->antena (independente do usuário)
    dsm = zeros(M, Lscat);
    ds1 = zeros(Lscat,1); % espalhador->antena de referência
    for l = 1:Lscat
        dxsl = URA(:,1) - scatterer_pos(l,1);
        dysl = URA(:,2) - scatterer_pos(l,2);
        dzsl = URA(:,3) - scatterer_pos(l,3);
        dsl  = sqrt(dxsl.^2 + dysl.^2 + dzsl.^2) + eps; % [M x 1]
        dsm(:,l) = dsl;
        ds1(l)   = dsl(1);
    end

    % Distâncias de 2 saltos (usuário->espalhador->antena m)
    dml = dsm + dus.';     % [M x Lscat], cada coluna l soma d_us(l) (broadcast)

    % beta_l_raw e campos raw por caminho l
    beta_l_raw = Calpha ./ (dml.^alpha);           % [M x Lscat] (potência)
    a_l_raw    = sqrt(beta_l_raw);                 % amplitudes
    % fase relativa (independente do usuário): (d_sm - d_sref)
    phi_l      = (2*pi/lambda) * (dsm - ds1.');    % [M x Lscat]
    h_l_raw    = a_l_raw .* exp(1j*phi_l) .* (ones(M,1)*Gamma.'); % [M x Lscat]

    % Potência "raw" agregada por caminho
    Pl_raw     = sum(abs(h_l_raw).^2, 1).';        % [Lscat x 1]
    Psum_raw   = sum(Pl_raw) + eps;

    % -------------------- Partição Rician (normalizações) ----------------
    % LoS normalizado: g0 = sqrt(K/(K+1)) (independe da posição!)
    g0   = sqrt( K_lin / (K_lin + 1) );
    h0   = g0 * h0_raw;

    % NLoS alvo: PN_target = 1/(K+1) * S_baseline
    PN_target = (1/(K_lin+1)) * S_baseline;

    % pesos w_l = Pl_raw / sum Pl_raw
    wl = Pl_raw ./ Psum_raw;                   % [Lscat x 1]
    target_l = PN_target .* wl;                % [Lscat x 1]
    gl = sqrt( target_l ./ (Pl_raw + eps) );   % [Lscat x 1]

    % Campo NLoS normalizado e total
    hN = h_l_raw * gl;                         % [M x 1], soma l (col) com pesos gl
    H  = h0 + hN;                              % canal final por elemento

    % -------------------- Derivadas dH/dx e dH/dy ------------------------
    % 1) Derivadas LoS (h0 = g0*h0_raw; g0 constante)
    % a0_raw = sqrt(C2)*dm^(-alpha/2);  phi0 = (2pi/lambda)*(dm - dref)
    uxm   = dxm ./ dm;                         % (x - x_m)/d_m
    uym   = dym ./ dm;

    da0dx = a0_raw .* ( -alpha/2 ) .* ( uxm ./ dm );      % [M x 1]
    da0dy = a0_raw .* ( -alpha/2 ) .* ( uym ./ dm );

    dphdx = (2*pi/lambda) * ( uxm - ( (x - x_ref)/dref ) );
    dphdy = (2*pi/lambda) * ( uym - ( (y - y_ref)/dref ) );

    dh0dx = g0 * ( exp(1j*phi0) .* ( da0dx + 1j .* a0_raw .* dphdx ) );
    dh0dy = g0 * ( exp(1j*phi0) .* ( da0dy + 1j .* a0_raw .* dphdy ) );

    % 2) Derivadas NLoS
    % hN = sum_l gl * h_l_raw(:,l)
    % ∂hN/∂x = sum_l [ (∂gl/∂x)*h_l_raw(:,l) + gl * ∂h_l_raw(:,l)/∂x ]
    %   - ∂phi_l/∂x = 0 (independe do usuário)
    %   - ∂a_l_raw/∂x via d_ml (∂d_ml/∂x = ∂d_us/∂x)

    % Derivadas de a_l_raw: a_l_raw = sqrt(Calpha) * (d_ml).^(-alpha/2) * Gamma(l)
    ddml_dx = repmat(ddus_dx.', M, 1);                     % [M x Lscat], mesmo por coluna
    ddml_dy = repmat(ddus_dy.', M, 1);

    da_l_dx = ( -alpha/2 ) .* ( a_l_raw ./ (dml + eps) ) .* ddml_dx;  % [M x Lscat]
    da_l_dy = ( -alpha/2 ) .* ( a_l_raw ./ (dml + eps) ) .* ddml_dy;

    % Derivadas h_l_raw (fase fixa):
    dh_l_raw_dx = da_l_dx .* exp(1j*phi_l) .* (ones(M,1)*Gamma.');    % [M x Lscat]
    dh_l_raw_dy = da_l_dy .* exp(1j*phi_l) .* (ones(M,1)*Gamma.');

    % Derivadas de gl (usar forma log para estabilidade)
    % gl = sqrt(target_l / Pl_raw) => ∂gl = 0.5*gl*( ∂log target_l - ∂log Pl_raw )
    % Derivadas de PN_target, Pl_raw e wl:
    % - dS_base/dx = sum_m d(beta0_raw)/dx  (pois S_base = sum beta0_raw)
    dbeta0dx = C2 * ( -alpha ) .* ( dm.^(-alpha-1) ) .* uxm;          % [M x 1]
    dbeta0dy = C2 * ( -alpha ) .* ( dm.^(-alpha-1) ) .* uym;

    dSbase_dx = sum(dbeta0dx);
    dSbase_dy = sum(dbeta0dy);

    dPN_dx = (1/(K_lin+1)) * dSbase_dx;
    dPN_dy = (1/(K_lin+1)) * dSbase_dy;

    % dPl_raw(l)/dx = |Gamma(l)|^2 * Calpha * sum_m (-alpha) d_ml^(-alpha-1) * ∂d_ml/∂x
    absG2 = abs(Gamma).^2;                                 % [Lscat x 1]
    dPl_dx = zeros(Lscat,1);
    dPl_dy = zeros(Lscat,1);
    for l = 1:Lscat
        t_dx = (-alpha) * sum( (dml(:,l)).^(-alpha-1) .* ddml_dx(:,l) );
        t_dy = (-alpha) * sum( (dml(:,l)).^(-alpha-1) .* ddml_dy(:,l) );
        dPl_dx(l) = absG2(l) * Calpha * t_dx;
        dPl_dy(l) = absG2(l) * Calpha * t_dy;
    end
    dPsum_dx = sum(dPl_dx) + eps;
    dPsum_dy = sum(dPl_dy) + eps;

    % dwl/dx, dwl/dy
    dw_dx = ( (dPl_dx .* Psum_raw) - (Pl_raw .* dPsum_dx) ) ./ (Psum_raw.^2);
    dw_dy = ( (dPl_dy .* Psum_raw) - (Pl_raw .* dPsum_dy) ) ./ (Psum_raw.^2);

    % d target_l / dθ
    dtarget_dx = dPN_dx .* wl + PN_target .* dw_dx;
    dtarget_dy = dPN_dy .* wl + PN_target .* dw_dy;

    % d log(target_l) / dθ e d log(Pl_raw) / dθ
    dlogT_dx = dtarget_dx ./ (target_l + eps);
    dlogT_dy = dtarget_dy ./ (target_l + eps);
    dlogP_dx = dPl_dx    ./ (Pl_raw   + eps);
    dlogP_dy = dPl_dy    ./ (Pl_raw   + eps);

    dgl_dx = 0.5 * gl .* ( dlogT_dx - dlogP_dx );         % [Lscat x 1]
    dgl_dy = 0.5 * gl .* ( dlogT_dy - dlogP_dy );

    % ∂hN/∂x = H_l_raw * dgl/dx  +  (∂H_l_raw/∂x) * gl
    dhN_dx = h_l_raw * dgl_dx + sum( dh_l_raw_dx .* (ones(M,1)*gl.'), 2 );
    dhN_dy = h_l_raw * dgl_dy + sum( dh_l_raw_dy .* (ones(M,1)*gl.'), 2 );

    % Gradientes finais
    dHdx = dh0dx + dhN_dx;     % [M x 1]
    dHdy = dh0dy + dhN_dy;

    % -------------------- FIM e CRB --------------------------------------
    Jxx = real( dHdx' * dHdx ) * pre_factor;
    Jyy = real( dHdy' * dHdy ) * pre_factor;
    Jxy = real( dHdx' * dHdy ) * pre_factor;

    FIM = [Jxx, Jxy; Jxy, Jyy];
    % Inversão estável
    detF = max(Jxx*Jyy - Jxy^2, eps);
    crb_x =  Jyy / detF;
    crb_y =  Jxx / detF;

    crb_euc = crb_x + crb_y;
end
