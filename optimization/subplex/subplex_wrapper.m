%[subplex_est, hist] = subplex_wrapper(URA, pos_est, Un, lambda, ref, x, y, 1e-5, 20);

function [subplex_est, history] = subplex_wrapper(URA, pos_init, Un, lambda, ref, x_lim, y_lim, tol, max_eval)

    % Função objetivo
    ps_func = @(xy) -music(xy(1), xy(2), URA, Un, lambda, ref);

    % Limites inferiores e superiores
    lb = [x_lim(1), y_lim(1)];
    ub = [x_lim(2), y_lim(2)];

    % Configurações do subplex
    options.maxeval = max_eval;  % número máximo de avaliações
    options.tol = tol;           % tolerância
    options.prt = 0;             % sem saída na tela
    options.bounds = [lb; ub];   % limites de busca
    options.restart = 0;         % sem restart automático

    % Chamada ao subplex (supondo sbplx.m disponível no path)
    [subplex_est, ~, history] = sbplx(ps_func, pos_init, options);
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
