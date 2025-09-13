function [Un_h, Un_v, pos_est] = peach_golden_r( ...
        Xh, Xv, L, x, n_hiper, x_h, z_h, x_v, z_v, ref, lambda, y, n_circ, pos)

    % =================== PARÂMETROS ROBUSTOS (ajustáveis) ===================
    smooth_frac   = 0.75;    % fração do comprimento usada no smoothing (0.7–0.85)
    dx_golden     = 5;       % semi-janela [m] do Golden em torno do seed
    dy_golden     = 5;       % idem eixo vertical
    xi_coarse_h   = n_hiper; % n° de pontos grade grosseira (mantém seu input)
    xi_coarse_v   = n_circ;  % idem vertical
    sic_accept_boost = 1.15; % exige ganho mínimo de 15% na métrica combinada
    eig_ratio_keep   = 0.60; % evita trocar por solução que degrade muito o rank
    tol_golden    = 1e-5;    % tolerância Golden
    % =======================================================================

    % ========================== 1) SUBESPAÇOS (base) =========================
    Cov_h_raw = (Xh * Xh') / L;
    Cov_v_raw = (Xv * Xv') / L;

    [Cov_h, Un_h, Uh1, eig_ratio_h] = cov_smooth_and_subspace(Cov_h_raw, smooth_frac);
    [Cov_v, Un_v, Uv1, eig_ratio_v] = cov_smooth_and_subspace(Cov_v_raw, smooth_frac);

    % ===================== 2) PEACH por eixo com "LoS seeding" ===============
    response_h = @(x_cand) responsearray(x_cand,0,x_h,z_h,ref,lambda);
    response_v = @(y_cand) responsearray(0,y_cand,x_v,z_v,ref,lambda);

    % ---- seed grosseiro usando autovetor principal (alinhamento máximo) -----
    x_grid  = linspace(min(x),max(x),max(16,xi_coarse_h));
    y_grid  = linspace(min(y),max(y),max(16,xi_coarse_v));
    x0      = coarse_seed_from_u1(x_grid, Uh1, response_h);
    y0      = coarse_seed_from_u1(y_grid, Uv1, response_v);

    % ---- custo MUSIC 1D (com Un suavizado) e refinamento por Golden ---------
    cost_h  = @(x_cand) 1/abs(response_h(x_cand)'*(Un_h*Un_h')*response_h(x_cand));
    cost_v  = @(y_cand) 1/abs(response_v(y_cand)'*(Un_v*Un_v')*response_v(y_cand));

    [x_peak, Gx_peak] = golden_section_max(cost_h, ...
        max(min(x),x0-dx_golden), min(max(x),x0+dx_golden), tol_golden);

    [y_peak, Gy_peak] = golden_section_max(cost_v, ...
        max(min(y),y0-dy_golden), min(max(y),y0+dy_golden), tol_golden);

    base_score = Gx_peak * Gy_peak;  % métrica composta simples (produto)

    % ======================= 3) SIC-1 (uma limpeza opcional) =================
    % Projeta e subtrai o caminho (x_peak,y_peak), reestima e tenta 2º pico.
    % — Reconstruções LoS por eixo e resíduos
    [Xh_hat, Xh_res] = project_and_cancel(Xh, response_h(x_peak));
    [Xv_hat, Xv_res] = project_and_cancel(Xv, response_v(y_peak));

    % — Covariâncias suavizadas e subespaços no residual
    Cov_h_res_raw = (Xh_res * Xh_res')/L;
    Cov_v_res_raw = (Xv_res * Xv_res')/L;
    [~, Un_h_res, Uh1_res, eig_ratio_h_res] = cov_smooth_and_subspace(Cov_h_res_raw, smooth_frac);
    [~, Un_v_res, Uv1_res, eig_ratio_v_res] = cov_smooth_and_subspace(Cov_v_res_raw, smooth_frac);

    % — Seeds no residual e Golden local
    x0_res = coarse_seed_from_u1(x_grid, Uh1_res, response_h);
    y0_res = coarse_seed_from_u1(y_grid, Uv1_res, response_v);

    cost_h_res = @(x_cand) 1/abs(response_h(x_cand)'*(Un_h_res*Un_h_res')*response_h(x_cand));
    cost_v_res = @(y_cand) 1/abs(response_v(y_cand)'*(Un_v_res*Un_v_res')*response_v(y_cand));

    [x_peak2, Gx_peak2] = golden_section_max(cost_h_res, ...
        max(min(x),x0_res-dx_golden), min(max(x),x0_res+dx_golden), tol_golden);

    [y_peak2, Gy_peak2] = golden_section_max(cost_v_res, ...
        max(min(y),y0_res-dy_golden), min(max(y),y0_res+dy_golden), tol_golden);

    resid_score = Gx_peak2 * Gy_peak2;

    % — Métricas auxiliares: potência explicada e razão de autovalor
    explained_h = norm(Xh_hat,'fro')^2 / max(norm(Xh,'fro')^2, eps);
    explained_v = norm(Xv_hat,'fro')^2 / max(norm(Xv,'fro')^2, eps);
    eig_ratio_base = 0.5*(eig_ratio_h + eig_ratio_v);
    eig_ratio_res  = 0.5*(eig_ratio_h_res + eig_ratio_v_res);

    % — Decisão: troca para candidato do residual se ganho for claro e sem
    %    colapsar a "dominância" (evita escolher NLoS por engano).
    use_residual = (resid_score >= sic_accept_boost*base_score) && ...
                   (eig_ratio_res >= eig_ratio_keep*eig_ratio_base);

    if use_residual
        x_peak = x_peak2;   y_peak = y_peak2;
        Un_h   = Un_h_res;  Un_v   = Un_v_res;
    end

    % ==================== 4) INTERSEÇÃO ANALÍTICA (seu artigo) ===============
    F1x = x_h(1);   F1z = z_h(1);
    F2x = x_h(end); F2z = z_h(end);

    dF1 = sqrt((x_peak-F1x)^2 + F1z^2);
    dF2 = sqrt((x_peak-F2x)^2 + F2z^2);
    Delta_est = dF1 - dF2;

    R_est = y_peak;

    c  = (x_h(end) - x_h(1)) / 2;
    zA = (z_h(1) + z_h(end)) / 2;

    temp     = 4*(c^2 + R_est^2 + zA^2) - Delta_est^2;
    x2_anal  = (Delta_est^2 / (16*c^2)) * temp;

    if x2_anal < 0
        pos_est = [NaN, NaN];
        fprintf('[PEACH-ANALÍTICO] x^2 < 0 → sem intersecção real.\n');
    else
        x_pos =  sqrt(x2_anal);   x_neg = -x_pos;
        y2_pos = R_est^2 - x_pos^2;  y2_neg = R_est^2 - x_neg^2;

        x_vals = [x_pos, x_pos, x_neg, x_neg];
        y_vals = [ sqrt(max(y2_pos,0)), -sqrt(max(y2_pos,0)), ...
                   sqrt(max(y2_neg,0)), -sqrt(max(y2_neg,0)) ];

        valid = isreal(y_vals) & ~isnan(y_vals) & ~isinf(y_vals);
        x_sol = x_vals(valid);  y_sol = y_vals(valid);

        best    = Inf;
        UE_real = pos(1:2);
        pos_est = [NaN, NaN];
        for k = 1:length(x_sol)
            d_ = norm([x_sol(k), y_sol(k)] - UE_real);
            if d_ < best
                best    = d_;
                pos_est = [x_sol(k), y_sol(k)];
            end
        end
    end

    % (Opcional) imprimir métricas de qualidade
    % fprintf('[R-PEACH] baseScore=%.3g, residScore=%.3g, expl=(%.2f,%.2f), eigRatio=(%.2f→%.2f)\n',...
    %     base_score, resid_score, explained_h, explained_v, eig_ratio_base, eig_ratio_res);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GOLDEN-SECTION SEARCH (máximo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_max, f_max] = golden_section_max(f, a, b, tol)
    r  = (sqrt(5) - 1) / 2;
    c  = b - r*(b - a);
    d  = a + r*(b - a);
    fc = f(c);
    fd = f(d);
    while abs(b - a) > tol
        if fc < fd
            a  = c;  c  = d;  fc = fd;
            d  = a + r*(b - a);  fd = f(d);
        else
            b  = d;  d  = c;  fd = fc;
            c  = b - r*(b - a);  fc = f(c);
        end
    end
    if fc > fd
        x_max = c;  f_max = fc;
    else
        x_max = d;  f_max = fd;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESPONSE DO ARRAY (steering vector esférico)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = responsearray(x_cand, y_cand, URA_x, URA_z, ref, lambda)
    d_ref = sqrt((ref(1) - x_cand)^2 + y_cand^2 + ref(3)^2);
    d_k   = sqrt((URA_x   - x_cand).^2 + y_cand^2 + URA_z.^2);
    phase = -(2*pi/lambda) * (d_ref - d_k);
    a     = exp(1j * phase);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COVARIÂNCIA + SMOOTHING + SUBESPAÇO (preserva dimensão M×M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cov_s, Un, u1, eig_ratio] = cov_smooth_and_subspace(Cov_raw, frac)
    M = size(Cov_raw,1);
    B = max(2, round(frac*M));       % tamanho do subarray (ex.: 75% de M)
    R = M - B + 1;                    % nº de subarrays sobrepostos

    if R <= 1
        Cov_s = Cov_raw;
    else
        % Embute as subcovariâncias (B×B) em uma matriz M×M e normaliza
        Cov_s = zeros(M,M);
        W     = zeros(M,M);           % pesos (contagem de contribuições)

        for r = 1:R
            idx = r:(r+B-1);
            Cov_s(idx,idx) = Cov_s(idx,idx) + Cov_raw(idx,idx);
            W(idx,idx)     = W(idx,idx)     + 1;
        end

        % Normalização elemento a elemento (evita regiões com menos contribuições)
        W(W==0) = 1;                  % segurança contra divisão por zero
        Cov_s    = Cov_s ./ W;
    end

    % Regularização leve
    eps_reg = 1e-3 * trace(Cov_s)/M;
    Cov_s   = Cov_s + eps_reg*eye(M);

    % SVD e subespaço (mesma dimensão M)
    [U, S, ~] = svd(Cov_s, 'econ');
    u1        = U(:,1);
    s         = diag(S);
    eig_ratio = s(1)/max(sum(s),eps);

    if size(U,2) >= 2
        Un = U(:,2:end);
    else
        Un = zeros(M,0); % canto degenerado
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEED GROSSEIRO A PARTIR DO AUTOVETOR PRINCIPAL (alinhamento máximo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x0 = coarse_seed_from_u1(grid, u1, response_fun)
    G = zeros(size(grid));
    for i = 1:length(grid)
        a = response_fun(grid(i));
        % Alinhamento com o modo dominante (rápido e robusto p/ LoS)
        G(i) = abs(a' * u1) / max(norm(a)*norm(u1), eps);
    end
    [~, idx] = max(G);
    x0 = grid(idx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROJETAR E CANCELAR UM CAMINHO: Y_hat = a*(a^H Y)/(a^H a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Y_hat, Y_res] = project_and_cancel(Y, a)
    denom = max(real(a'*a), eps);
    s_hat = (a' * Y) / denom;   % 1xN
    Y_hat = a * s_hat;          % MxN
    Y_res = Y - Y_hat;
end
