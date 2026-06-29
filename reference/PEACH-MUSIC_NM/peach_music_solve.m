function [Un_h, Un_v, pos_est, best_dist] = peach_music_solve(Xh, Xv, L, x, n_hiper, ...
    x_h, z_h, x_v, z_v, x_ref_h, z_ref_h, x_ref_v, z_ref_v, ...
    lambda, y, n_circ, pos)

    %-----------------------------------------------
    % DIVISÃO DO SUBESPAÇO (Vn)
    %-----------------------------------------------
    Cov_h = (Xh * Xh') / L;
    Cov_v = (Xv * Xv') / L;

    [Uh, ~, ~] = svd(Cov_h);
    [Uv, ~, ~] = svd(Cov_v);
    
    Un_h = Uh(:, 2:end);
    Un_v = Uv(:, 2:end);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Subarranjo horizontal => HIPÉRBOLE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    response_h = @(x_cand) responsearray(x_cand, 0, x_h, z_h, x_ref_h, z_ref_h, lambda);
    x_candidates = linspace(-max(x), max(x), n_hiper);
    [~, x_peak] = music(x_candidates, Un_h, response_h);
    
    F1x = x_h(1);  F1z = z_h(1);
    F2x = x_h(end); F2z = z_h(end);

    dF1 = sqrt((x_peak - F1x)^2 + F1z^2);
    dF2 = sqrt((x_peak - F2x)^2 + F2z^2);
    Delta_est = dF1 - dF2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Subarranjo vertical => CÍRCULO
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    response_v = @(y_cand) responsearray(0, y_cand, x_v, z_v, x_ref_v, z_ref_v, lambda);
    y_candidates = linspace(0, max(y), n_circ);
    [~, y_peak] = music(y_candidates, Un_v, response_v);
    R_est = y_peak;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INTERSEÇÃO: HIPÉRBOLE vs CÍRCULO (via solve)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    syms xs ys real
    expr_d1 = sqrt((xs - F1x)^2 + ys^2 + F1z^2);
    expr_d2 = sqrt((xs - F2x)^2 + ys^2 + F2z^2);
    eq_hip  = expr_d1 - expr_d2 - Delta_est;
    eq_circ = xs^2 + ys^2 - R_est^2;

    sol_ = solve([eq_hip, eq_circ], [xs, ys], 'Real', true);
    x_sol_ = double(sol_.xs);
    y_sol_ = double(sol_.ys);

    if isempty(x_sol_)
        pos_est = [NaN, NaN];
        best_dist = NaN;
        fprintf('[PEACH-MUSIC] Sem soluções reais.\n');
    else
        best_dist = Inf;
        UE_real = pos(1:2);
        for kk = 1:length(x_sol_)
            cand = [x_sol_(kk), y_sol_(kk)];
            d_ = norm(cand - UE_real);
            if d_ < best_dist
                best_dist = d_;
                pos_est = cand;
            end
        end
    end
    fprintf('PEACH final = (%.3f, %.3f)\n', pos_est(1), pos_est(2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNÇÃO LOCAL - RESPONSE ARRAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function a = responsearray(x_cand, y_cand, URA_x, URA_z, x_ref, z_ref, lambda)
    M = numel(URA_x);
    a = zeros(M,1);

    d_ref = sqrt((x_ref - x_cand)^2 + y_cand^2 + z_ref^2);
    for m = 1:M
        x_k = URA_x(m);
        z_k = URA_z(m);
        d_k = sqrt((x_k - x_cand)^2 + y_cand^2 + z_k^2);
        phase = -(2*pi/lambda)*(d_ref - d_k);
        a(m) = exp(1j * phase);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNÇÃO LOCAL - MUSIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [peak_val, peak_coord] = music(candidates, Un, response_fun)
    G = zeros(size(candidates));
    for i = 1:length(candidates)
        a = response_fun(candidates(i));
        G(i) = 1 / abs(a' * (Un * Un') * a);
    end
    [peak_val, idx_max] = max(G);
    peak_coord = candidates(idx_max);
end