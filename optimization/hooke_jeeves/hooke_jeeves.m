% % Hooke–Jeeves refinement
% [ hj_est, hj_history ] = hooke_jeeves(URA, est_peach, Un, lambda, ref, ...
%     deltaArea, numIterNM, tol, x, y, true);

% hooke_jeeves.m
function [hj_est, history] = hooke_jeeves(URA, x0, Un, lambda, ref, delta0, max_iter, tol, ...
    x_bounds, y_bounds, verbose)

    % Derivative-free pattern search (Hooke–Jeeves)
    obj = @(xy) -music(xy(1), xy(2), URA, Un, lambda, ref);  % minimiza -PMUSIC
    x_curr = x0(:)';
    f_curr = obj(x_curr);
    delta  = delta0;
    history = x_curr;
    iter = 0;

    while delta > tol && iter < max_iter
        % 1) Exploração coordenada
        x_exp = x_curr;
        for i = 1:2
            % avanço
            xt = x_exp;
            xt(i) = min(max(x_exp(i) + delta, x_bounds(1)), x_bounds(2));
            if obj(xt) < obj(x_exp)
                x_exp = xt;
            else
                % recuo
                xt = x_exp;
                xt(i) = min(max(x_exp(i) - delta, x_bounds(1)), x_bounds(2));
                if obj(xt) < obj(x_exp)
                    x_exp = xt;
                end
            end
        end

        % 2) Pattern move ou redução de passo
        if obj(x_exp) < f_curr
            x_new = 2*x_exp - x_curr;
            % limites
            x_new(1) = min(max(x_new(1), x_bounds(1)), x_bounds(2));
            x_new(2) = min(max(x_new(2), y_bounds(1)), y_bounds(2));
            x_curr = x_new;
            f_curr = obj(x_curr);
        else
            delta = delta / 2;
        end

        iter = iter + 1;
        history(end+1, :) = x_curr;  %#ok<AGROW>
    end

    hj_est = x_curr;
    if verbose
        %fprintf('Hooke–Jeeves convergiu em %d iterações, Δ = %.3e\n', iter, delta);
    end
end

% subfunção local: cálculo do pseudospectro MUSIC
function valPS = music(x, y, URA, Un, lambda, ref)
    d_ref = sqrt((ref(1)-x)^2 + (ref(2)-y)^2 + ref(3)^2);
    d_km  = sqrt((URA(:,1)-x).^2 + (URA(:,2)-y).^2 + URA(:,3).^2);
    phase = -(2*pi/lambda) * (d_ref - d_km);
    a     = exp(1j * phase);
    den   = abs(a' * (Un*Un') * a);
    % substitui o '? :' por if/else
    if den < 1e-12
        valPS = 1e6;
    else
        valPS = 1/den;
    end
end