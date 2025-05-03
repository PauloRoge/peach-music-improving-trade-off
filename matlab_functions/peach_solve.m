function [Un_h, Un_v, pos_est, best_dist] = peach_solve(Xh, Xv, L, x, n_hiper, URA, ref, lambda, y, n_circ, pos)
    %------------------------------------------------
    % 1) Subespaço ruído para subarranjos H e V
    Cov_h = (Xh * Xh') / L;
    Cov_v = (Xv * Xv') / L;
    [Uh, ~, ~] = svd(Cov_h);
    [Uv, ~, ~] = svd(Cov_v);
    Un_h = Uh(:, 2:end);
    Un_v = Uv(:, 2:end);

    %------------------------------------------------
    % 2) Extrai subarranjo horizontal (mesma altura z_ref)
    hor_idx = abs(URA(:,3) - ref(3)) < 1e-8;
    x_h = URA(hor_idx, 1);
    z_h = URA(hor_idx, 3);

    % Subarranjo vertical (mesma coordenada x_ref)
    ver_idx = abs(URA(:,1) - ref(1)) < 1e-8;
    x_v = URA(ver_idx, 1);
    z_v = URA(ver_idx, 3);

    %------------------------------------------------
    % 3) Estima Δ pela hipérbole (x_peak)
    response_h = @(x_cand) responsearray(x_cand, 0, x_h, z_h, ref(1), ref(3), lambda);
    x_candidates = linspace(-max(x), max(x), n_hiper);
    [~, x_peak] = music(x_candidates, Un_h, response_h);

    F1x = x_h(1);    F1z = z_h(1);
    F2x = x_h(end);  F2z = z_h(end);
    dF1 = sqrt((x_peak - F1x)^2 + F1z^2);
    dF2 = sqrt((x_peak - F2x)^2 + F2z^2);
    Delta_est = dF1 - dF2;

    %------------------------------------------------
    % 4) Estima R pelo círculo (y_peak)
    response_v = @(y_cand) responsearray(0, y_cand, x_v, z_v, ref(1), ref(3), lambda);
    y_candidates = linspace(0, max(y), n_circ);
    [~, y_peak] = music(y_candidates, Un_v, response_v);
    R_est = y_peak;

    %------------------------------------------------
    % 5) Interseção: resolve hipérbole vs círculo
    syms xs ys real
    expr_d1 = sqrt((xs - F1x)^2 + ys^2 + F1z^2);
    expr_d2 = sqrt((xs - F2x)^2 + ys^2 + F2z^2);
    sol = solve([expr_d1 - expr_d2 - Delta_est, xs^2 + ys^2 - R_est^2], [xs, ys], 'Real', true);
    x_sol = double(sol.xs);  y_sol = double(sol.ys);

    if isempty(x_sol)
        pos_est   = [NaN, NaN];
        best_dist = NaN;
        %fprintf('[PEACH-MUSIC] Sem soluções reais.\n');
    else
        best_dist = Inf;
        UE = pos(1:2);
        for k = 1:length(x_sol)
            cand = [x_sol(k), y_sol(k)];
            d = norm(cand - UE);
            if d < best_dist
                best_dist = d;
                pos_est   = cand;
            end
        end
    end
    %fprintf('PEACH final = (%.3f, %.3f)\n', pos_est(1), pos_est(2));
end

% ================================================
% FUNÇÃO LOCAL: calcula steering vector em 2D
function a = responsearray(x_cand, y_cand, URA_x, URA_z, x_ref, z_ref, lambda)
    M = numel(URA_x);
    a = zeros(M,1);
    d_ref = sqrt((x_ref - x_cand)^2 + y_cand^2 + z_ref^2);
    for m = 1:M
        d_k = sqrt((URA_x(m) - x_cand)^2 + y_cand^2 + URA_z(m)^2);
        a(m) = exp(-1j*(2*pi/lambda)*(d_ref - d_k));
    end
end

% ================================================
% FUNÇÃO LOCAL: encontra pico do pseudo-espectro
function [peak_val, peak_coord] = music(cands, Un, resp_fun)
    G = arrayfun(@(c) 1/abs(resp_fun(c)' * (Un*Un') * resp_fun(c)), cands);
    [peak_val, idx] = max(G);
    peak_coord = cands(idx);
end
