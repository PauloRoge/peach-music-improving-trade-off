% [subplex_est, history, total_evals] = subplex_wrapper_test(URA, pos_init, Un, ...
%     lambda, ref, x, y, tol, 560);

function [subplex_est, history, total_evals] = subplex_wrapper_test(URA, pos_init, Un, ...
    lambda, ref, x_lim, y_lim, tol, max_eval)
% SUBPLEX_WRAPPER  Invoca sbplx e retorna o total de avaliações
%   [subplex_est, history, total_evals] = subplex_wrapper(...)
%   total_evals: número total de vezes que a função objetivo foi avaliada.

    % Função objetivo (negação do MUSIC para maximização)
    ps_func = @(xy) -music(xy(1), xy(2), URA, Un, lambda, ref);

    % Ajusta limites de busca
    lb = [x_lim(1), y_lim(1)];
    ub = [x_lim(2), y_lim(2)];

    % Configurações do Subplex
    options.maxeval = max_eval;     % avaliação máxima permitida
    options.tol     = tol;          % tolerância de convergência
    options.prt     = 0;            % sem saída textual
    options.bounds  = [lb; ub];     % limites inferiores e superiores
    options.restart = 0;            % sem restart automático

    % Chama sbplx e captura neval
    [subplex_est, ~, history, neval] = sbplx_test(ps_func, pos_init, options);

    % Aqui somaríamos múltiplos neval caso houvesse várias chamadas,
    % mas como sbplx é chamado apenas uma vez, total_evals = neval.
    total_evals = neval;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCAO LOCAL - MUSIC PSEUDOESPECTRO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function valPS = music(x, y, URA, Un, lambda, ref)
    d_ref = sqrt((ref(1) - x)^2 + (ref(2) - y)^2 + ref(3)^2);
    d_km = sqrt((URA(:,1) - x).^2 + (URA(:,2) - y).^2 + URA(:,3).^2);
    phase_diff = -(2*pi/lambda) * (d_ref - d_km);
    a = exp(1j * phase_diff);
    den = abs(a' * (Un * Un') * a);
    if den < 1e-12
        valPS = 1e6;
    else
        valPS = 1 / den;
    end
end