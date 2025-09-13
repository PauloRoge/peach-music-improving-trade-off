function [Un_h, Un_v, pos_est] = peach_rand(Xh, Xv, L, x, n_hiper, ...
    x_h, z_h, x_v, z_v, ref, ...
    lambda, y, n_circ, pos)

    %-----------------------------------------------
    % SUBESPACOS COM DIAGONAL LOADING
    %-----------------------------------------------
    Cov_h = (Xh * Xh') / L;                    % covariância estimada
    Cov_v = (Xv * Xv') / L;
    
    % escolhe ε como uma fração (e.g. 10^{-3}) da potência média
    M_h = size(Cov_h,1);
    M_v = size(Cov_v,1);
    eps_h = 1e-3 * trace(Cov_h) / M_h;
    eps_v = 1e-3 * trace(Cov_v) / M_v;
    
    Cov_h = Cov_h + eps_h * eye(M_h);          % aplicação do diagonal loading
    Cov_v = Cov_v + eps_v * eye(M_v);
    
    [Uh, ~, ~] = svd(Cov_h);                    % SVD regularizada
    [Uv, ~, ~] = svd(Cov_v);
    
    Un_h = Uh(:, 2:end);
    Un_v = Uv(:, 2:end);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Subarranjo horizontal => HIPERBOLE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    response_h = @(x_cand) responsearray(x_cand, 0, x_h, z_h, ref, lambda);
    %K = 96;  % número de pontos iniciais
    x_candidates = min(x) + (max(x) - min(x)) * rand(1, n_hiper);
    [~, x_peak] = music(x_candidates, Un_h, response_h);

    F1x = x_h(1);  F1z = z_h(1);
    F2x = x_h(end); F2z = z_h(end);

    dF1 = sqrt((x_peak - F1x)^2 + F1z^2);
    dF2 = sqrt((x_peak - F2x)^2 + F2z^2);
    Delta_est = dF1 - dF2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Subarranjo vertical => CIRCULO
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    response_v = @(y_cand) responsearray(0, y_cand, x_v, z_v, ref, lambda);
    %K = 48;
    y_candidates = min(y) + (max(y) - min(y)) * rand(1, n_circ);
    % y_candidates = linspace(0, max(y), n_circ);
    [~, y_peak] = music(y_candidates, Un_v, response_v);
    R_est = y_peak;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INTERSECAO ANALITICA (conforme artigo)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    c = (F2x - F1x) / 2;
    zA = (F1z + F2z) / 2;  % media das alturas dos focos (se forem diferentes)

    if c <= 0
        error('Parametro c invalido. Esperado F2x > F1x.');
    end

    % Calculo de x^2 (analitico) com base na deducao do artigo
    temp = 4 * (c^2 + R_est^2 + zA^2) - Delta_est^2;
    x2_anal = (Delta_est^2 / (16 * c^2)) * temp;

    if x2_anal < 0
        x_sol = []; y_sol = [];
        fprintf('[PEACH-MUSIC Analitico] x^2 < 0 => sem intersecao real.\n');
    else
        x_pos = sqrt(x2_anal);
        x_neg = -x_pos;

        y2_pos = R_est^2 - x_pos^2;
        y2_neg = R_est^2 - x_neg^2;

        y_vals = [sqrt(max(y2_pos, 0)), -sqrt(max(y2_pos, 0)), ...
                  sqrt(max(y2_neg, 0)), -sqrt(max(y2_neg, 0))];
        x_vals = [x_pos, x_pos, x_neg, x_neg];

        valid_idx = isreal(y_vals) & ~isnan(y_vals) & ~isinf(y_vals);
        x_sol = x_vals(valid_idx);
        y_sol = y_vals(valid_idx);
    end

    % Selecao da melhor solucao
    if isempty(x_sol)
         pos_est = [NaN, NaN];
        % best_dist = NaN;
        fprintf('[PEACH-MUSIC] Sem solucoes reais.\n');
    else
        best_dist = Inf;
        UE_real = pos(1:2);
        for kk = 1:length(x_sol)
            cand = [x_sol(kk), y_sol(kk)];
            d_ = norm(cand - UE_real);
            if d_ < best_dist
                best_dist = d_;
                pos_est = cand;
            end
        end
    end
    % fprintf('PEACH final = (%.3f, %.3f)\n', pos_est(1), pos_est(2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCAO LOCAL - RESPONSE ARRAY CANDIDATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function a = responsearray(x_cand, y_cand, URA_x, URA_z, ref, lambda)
    % distância referência–usuário
    d_ref = sqrt((ref(1) - x_cand).^2 + y_cand.^2 + ref(3).^2);

    % distâncias de todos os elementos de uma só vez
    d_k   = sqrt((URA_x - x_cand).^2 + y_cand.^2 + URA_z.^2);

    % fase vetorizada
    phase = -(2*pi/lambda) * (d_ref - d_k);

    % steering vector
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
