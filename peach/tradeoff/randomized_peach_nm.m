function [Un_h, Un_v, pos_est] = randomized_peach_nm(Xh, Xv, L, x, ...
    x_h, z_h, x_v, z_v, ref, ...
    lambda, y, pos)

    %-----------------------------------------------
    % SUBESPACOS COM DIAGONAL LOADING
    %-----------------------------------------------
    Cov_h = (Xh * Xh') / L;
    Cov_v = (Xv * Xv') / L;
    M_h   = size(Cov_h,1);
    M_v   = size(Cov_v,1);
    eps_h = 1e-3 * trace(Cov_h) / M_h;
    eps_v = 1e-3 * trace(Cov_v) / M_v;
    Cov_h = Cov_h + eps_h * eye(M_h);
    Cov_v = Cov_v + eps_v * eye(M_v);
    [Uh, ~, ~] = svd(Cov_h);
    [Uv, ~, ~] = svd(Cov_v);
    Un_h = Uh(:, 2:end);
    Un_v = Uv(:, 2:end);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Subarranjo horizontal => HIPERBOLE via randomized + Nelder – Mead
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % define subarranjo horizontal
    URA_h = [ x_h(:), zeros(numel(x_h),1), z_h(:) ];
    response_h = @(xc) responsearray(xc, 0, x_h, z_h, ref, lambda);

    % 1) amostragem aleatória
    K       = 96;
    x_rand  = min(x) + (max(x)-min(x)) * rand(K,1);
    G_rand  = zeros(K,1);
    for k = 1:K
        a         = response_h(x_rand(k));
        G_rand(k) = 1/abs(a' * (Un_h*Un_h') * a);
    end

    % 2) melhor semente
    [~, idxb1] = max(G_rand);
    seed_best = x_rand(idxb1);

    % 3) refinamento Nelder–Mead
    delta_x     = (max(x)-min(x)) / 20;
    max_iter_nm = 50;
    tol_nm      = 1e-6;
    x_bounds    = [min(x), max(x)];
    y_bounds_h  = [0, 0];
    
    [nm_est_h, ~] = nelder_mead( ...
        URA_h, [seed_best, 0], Un_h, lambda, ref, ...
        delta_x, max_iter_nm, tol_nm, false, ...
        x_bounds, y_bounds_h);
    x_peak = nm_est_h(1);

    % calcula Δ_est
    F1x      = x_h(1);   F1z = z_h(1);
    F2x      = x_h(end); F2z = z_h(end);
    dF1      = sqrt((x_peak - F1x)^2 + F1z^2);
    dF2      = sqrt((x_peak - F2x)^2 + F2z^2);
    Delta_est = dF1 - dF2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Subarranjo vertical => CÍRCULO via randomized + Nelder – Mead
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    URA_v     = [ x_v(:), zeros(numel(x_v),1), z_v(:) ];
    response_v = @(yc) responsearray(0, yc, x_v, z_v, ref, lambda);

    % 1) amostragem aleatória
    K2      = 48;
    y_rand  = min(y) + (max(y)-min(y)) * rand(K2,1);
    G2      = zeros(K2,1);
    for k = 1:K2
        a2      = response_v(y_rand(k));
        G2(k)   = 1/abs(a2' * (Un_v*Un_v') * a2);
    end

    % 2) melhor semente
    [~, idxb2] = max(G2);
    seed2_best = y_rand(idxb2);

    % 3) refinamento Nelder–Mead
    delta_y    = (max(y)-min(y)) / 20;
    x_bounds_v = [0, 0];
    y_bounds   = [min(y), max(y)];
    [nm_est_v, ~] = nelder_mead( ...
        URA_v, [0, seed2_best], Un_v, lambda, ref, ...
        delta_y, max_iter_nm, tol_nm, false, ...
        x_bounds_v, y_bounds);
    R_est = nm_est_v(2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INTERSEÇÃO ANALÍTICA (PEACH-MUSIC)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    c       = (F2x - F1x) / 2;
    zA      = (F1z + F2z) / 2;
    temp    = 4*(c^2 + R_est^2 + zA^2) - Delta_est^2;
    x2_anal = (Delta_est^2/(16*c^2)) * temp;
    if x2_anal < 0
        pos_est = [NaN, NaN];
    else
        x_sol = [ sqrt(x2_anal), -sqrt(x2_anal) ];
        y_sol = sqrt(max(R_est^2 - x_sol.^2, 0));
        UE    = pos(1:2);
        cands = [ x_sol(:), y_sol(:) ];
        dists = vecnorm(cands - UE, 2, 2);
        [~, idx] = min(dists);
        pos_est  = cands(idx, :);
    end
end

%% função local: steering near-field
function a = responsearray(x_cand, y_cand, URA_x, URA_z, ref, lambda)
    d_ref = sqrt((ref(1)-x_cand)^2 + y_cand^2 + ref(3)^2);
    d_k   = sqrt((URA_x - x_cand).^2 + y_cand.^2 + URA_z.^2);
    phase = -(2*pi/lambda) * (d_ref - d_k);
    a     = exp(1j * phase);
end
