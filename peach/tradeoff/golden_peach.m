function [Un_h, Un_v, pos_est] = golden_peach(Xh, Xv, L, x, ...
    x_h, z_h, x_v, z_v, ref, ...
    lambda, y, pos)

    %-----------------------------------------------
    % 1) SUBESPACOS COM DIAGONAL LOADING
    %-----------------------------------------------
    Cov_h = (Xh * Xh') / L;
    Cov_v = (Xv * Xv') / L;
    
    M_h = size(Cov_h,1); 
    M_v = size(Cov_v,1);
    eps_h = 1e-3 * trace(Cov_h) / M_h;
    eps_v = 1e-3 * trace(Cov_v) / M_v;
    
    Cov_h = Cov_h + eps_h * eye(M_h);
    Cov_v = Cov_v + eps_v * eye(M_v);
    
    [Uh,~,~] = svd(Cov_h);
    [Uv,~,~] = svd(Cov_v);
    
    Un_h = Uh(:,2:end);
    Un_v = Uv(:,2:end);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2) SUBARRANJO HORIZONTAL ⇒ HIPÉRBOLE (Golden‐section Search)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    response_h = @(x_cand) responsearray(x_cand, 0, x_h, z_h, ref, lambda);
    cost_h     = @(x_cand) 1/abs( response_h(x_cand)' * (Un_h*Un_h') * response_h(x_cand) );
    [x_peak,~] = golden_section_max(cost_h, min(x), max(x), 1e-3);

    F1x = x_h(1);   F1z = z_h(1);
    F2x = x_h(end); F2z = z_h(end);

    dF1 = sqrt((x_peak - F1x)^2 + F1z^2);
    dF2 = sqrt((x_peak - F2x)^2 + F2z^2);
    Delta_est = dF1 - dF2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3) SUBARRANJO VERTICAL ⇒ CÍRCULO (Golden‐section Search)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    response_v = @(y_cand) responsearray(0, y_cand, x_v, z_v, ref, lambda);
    cost_v     = @(y_cand) 1/abs( response_v(y_cand)' * (Un_v*Un_v') * response_v(y_cand) );
    [y_peak,~] = golden_section_max(cost_v, min(y), max(y), 1e-3);
    R_est      = y_peak;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4) INTERSEÇÃO ANALÍTICA (conforme artigo)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    c  = (x_h(end) - x_h(1)) / 2;
    zA = (z_h(1) + z_h(end)) / 2;

    temp     = 4*(c^2 + R_est^2 + zA^2) - Delta_est^2;
    x2_anal  = (Delta_est^2 / (16*c^2)) * temp;

    if x2_anal < 0
        pos_est = [NaN, NaN];
        fprintf('[PEACH-ANALÍTICO] x^2 < 0 → sem intersecção real.\n');
    else
        x_pos =  sqrt(x2_anal);
        x_neg = -x_pos;

        y2_pos = R_est^2 - x_pos^2;
        y2_neg = R_est^2 - x_neg^2;

        x_vals = [x_pos, x_pos, x_neg, x_neg];
        y_vals = [ sqrt(max(y2_pos,0)), -sqrt(max(y2_pos,0)), ...
                   sqrt(max(y2_neg,0)), -sqrt(max(y2_neg,0)) ];

        valid = isreal(y_vals) & ~isnan(y_vals) & ~isinf(y_vals);
        x_sol = x_vals(valid);
        y_sol = y_vals(valid);

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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNÇÃO GOLDEN‐SECTION SEARCH PARA MÁXIMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_max, f_max] = golden_section_max(f, a, b, tol)
    r  = (sqrt(5) - 1) / 2;
    c  = b - r*(b - a);
    d  = a + r*(b - a);
    fc = f(c);
    fd = f(d);

    while abs(b - a) > tol
        if fc < fd
            a  = c;
            c  = d;
            fc = fd;
            d  = a + r*(b - a);
            fd = f(d);
        else
            b  = d;
            d  = c;
            fd = fc;
            c  = b - r*(b - a);
            fc = f(c);
        end
    end

    if fc > fd
        x_max = c;
        f_max = fc;
    else
        x_max = d;
        f_max = fd;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNÇÃO DE RESPONSE DO ARRAY (roda steering‐vector)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = responsearray(x_cand, y_cand, URA_x, URA_z, ref, lambda)
    d_ref = sqrt((ref(1) - x_cand)^2 + y_cand^2 + ref(3)^2);
    d_k   = sqrt((URA_x   - x_cand).^2 + y_cand^2 + URA_z.^2);
    phase = -(2*pi/lambda) * (d_ref - d_k);
    a     = exp(1j * phase);
end
