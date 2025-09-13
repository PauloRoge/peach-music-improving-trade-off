function [Un_h, Un_v, pos_est] = randomized_peach(Yh, Yv, L, x, n_hiper, ...
    x_h, z_h, x_v, z_v, ref, ...
    lambda, y, n_circ, pos)

    %-----------------------------------------------
    % SUBESPACOS COM DIAGONAL LOADING
    %-----------------------------------------------
    Cov_h = (Yh * Yh') / L;                    % covariância estimada
    Cov_v = (Yv * Yv') / L;
    
    % escolhe ε como uma fração (e.g. 10^{-3}) da potência média
    M_h = size(Cov_h,1);
    M_v = size(Cov_v,1);
    eps_h = 1e-3 * trace(Cov_h) / M_h;
    eps_v = 1e-3 * trace(Cov_v) / M_v;
    
    % aplicação do diagonal loading
    Cov_h = Cov_h + eps_h * eye(M_h);
    Cov_v = Cov_v + eps_v * eye(M_v);

    % fatoração SVD regularizada
    [Uh, ~, ~] = svd(Cov_h);
    [Uv, ~, ~] = svd(Cov_v);
    
    % subespaços de ruído
    Un_h = Uh(:, 2:end);
    Un_v = Uv(:, 2:end);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Subarranjo horizontal => HIPERBOLE via randomized-MUSIC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    response_h = @(xc) responsearray(xc, 0, x_h, z_h, ref, lambda);
    % 1) amostragem aleatória
    K = n_hiper;  % número de pontos iniciais
    x_rand = min(x) + (max(x) - min(x)) * rand(1, K);
    G_rand = zeros(K,1);
    for k = 1:K
        a = response_h(x_rand(k));
        G_rand(k) = 1/abs(a' * (Un_h*Un_h') * a);
    end
    % 2) seleciona P sementes
    P = 10;
    [~, idxs] = maxk(G_rand, P);
    seeds = x_rand(idxs);
    % 3) refinamento coarse-to-fine
    x_ref = zeros(P,1);
    dx = (max(x)-min(x))/20;
    for p = 1:P
        xi = linspace(seeds(p)-dx, seeds(p)+dx, round(n_hiper/5));
        xi = max(min(xi, max(x)), min(x));
        [~, best] = music(xi, Un_h, response_h);
        x_ref(p) = best;
    end
    % 4) escolhe o melhor refinado
    G_ref = zeros(P,1);
    for p = 1:P
        a = response_h(x_ref(p));
        G_ref(p) = 1/abs(a' * (Un_h*Un_h') * a);
    end
    [~, ibest] = max(G_ref);
    x_peak = x_ref(ibest);

    % cálculos de Δ_est
    F1x = x_h(1);    F1z = z_h(1);
    F2x = x_h(end);  F2z = z_h(end);
    dF1 = sqrt((x_peak - F1x)^2 + F1z^2);
    dF2 = sqrt((x_peak - F2x)^2 + F2z^2);
    Delta_est = dF1 - dF2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Subarranjo vertical => CÍRCULO via randomized-MUSIC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    response_v = @(yc) responsearray(0, yc, x_v, z_v, ref, lambda);
    % 1) amostragem aleatória
    K2 = n_circ;
    y_rand = min(y) + (max(y) - min(y)) * rand(1, K2);
    G2 = zeros(K2,1);
    for k = 1:K2
        a = response_v(y_rand(k));
        G2(k) = 1/abs(a' * (Un_v*Un_v') * a);
    end
    % 2) seleciona P2 sementes
    P2 = 10;
    [~, id2] = maxk(G2, P2);
    seeds2 = y_rand(id2);
    % 3) refinamento coarse-to-fine
    y_ref = zeros(P2,1);
    dy = (max(y)-min(y))/20;
    for p = 1:P2
        yi = linspace(seeds2(p)-dy, seeds2(p)+dy, round(n_circ/5));
        yi = max(min(yi, max(y)), min(y));
        [~, best] = music(yi, Un_v, response_v);
        y_ref(p) = best;
    end
    % 4) escolhe o melhor refinado
    G2r = zeros(P2,1);
    for p = 1:P2
        a = response_v(y_ref(p));
        G2r(p) = 1/abs(a' * (Un_v*Un_v') * a);
    end
    [~, ib2] = max(G2r);
    R_est = y_ref(ib2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INTERSECAO ANALITICA (conforme artigo)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    c  = (F2x - F1x) / 2;
    zA = (F1z + F2z) / 2;  % média das alturas dos focos
    if c <= 0
        error('Parâmetro c inválido. Esperado F2x > F1x.');
    end
    temp   = 4 * (c^2 + R_est^2 + zA^2) - Delta_est^2;
    x2_anal = (Delta_est^2 / (16 * c^2)) * temp;
    if x2_anal < 0
        pos_est = [NaN, NaN];
    else
        x_sol = [ sqrt(x2_anal), -sqrt(x2_anal) ];
        y_sol = sqrt(max(R_est^2 - x_sol.^2, 0));
        
        % escolhe a solução mais próxima da posição real
        UE    = pos(1:2);
        % monta corretamente cada candidato como linha [x y]
        cands = [ x_sol(:), y_sol(:) ];  
        % distância euclidiana de cada linha [x y] até UE
        dists = vecnorm(cands - UE, 2, 2);
        [~, idx] = min(dists);
        pos_est = cands(idx, :);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCAO LOCAL - RESPONSE ARRAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function a = responsearray(x_cand, y_cand, URA_x, URA_z, ref, lambda)
    d_ref = sqrt((ref(1) - x_cand)^2 + y_cand^2 + ref(3)^2);
    d_k   = sqrt((URA_x - x_cand).^2 + y_cand.^2 + URA_z.^2);
    phase = -(2*pi/lambda) * (d_ref - d_k);
    a     = exp(1j * phase);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCAO LOCAL - MUSIC CANDIDATES
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
