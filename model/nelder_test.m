function [nm_est, simplex_history] = nelder_test(URA, pos_peach, Un, lambda, ref, ...
                                                delta, max_iter, tol, verbose, x, y, ...
                                                alpha, gamma, rho, sigma)

    if nargin < 8,  tol     = 1e-6;  end
    if nargin < 9,  verbose = false; end
    if nargin < 12, alpha = 1.0; end
    if nargin < 13, gamma = 2.0; end
    if nargin < 14, rho   = 0.5; end
    if nargin < 15, sigma = 0.5; end

    x_min = min(x(:)); x_max = max(x(:));
    y_min = min(y(:)); y_max = max(y(:));

    % função objetivo
    ps_func = @(x,y) music(x, y, URA, Un, lambda, ref);
    nm_val  = @(xy) -ps_func(xy(1), xy(2));

    % simplex inicial
    p1 = pos_peach;
    p2 = p1 + [delta, 0];
    p3 = p1 + [0, delta];

    % restrição ao espaço de busca
    p2(1) = min(max(p2(1), x_min), x_max);
    p2(2) = min(max(p2(2), y_min), y_max);
    p3(1) = min(max(p3(1), x_min), x_max);
    p3(2) = min(max(p3(2), y_min), y_max);

    simplexNM = [p1; p2; p3];
    fvals = zeros(3,1);

    for si = 1:3
        fvals(si) = nm_val(simplexNM(si,:));
    end

    simplex_history = simplexNM;

    for iterNM = 1:max_iter
        % ordena pelo valor da função
        [fvals, idx] = sort(fvals);
        simplexNM = simplexNM(idx,:);

        % checa convergência
        f_diff = max(abs(fvals - fvals(1)));
        if f_diff < tol
            break;
        end

        % centróide dos melhores dois
        mid_ = (simplexNM(1,:) + simplexNM(2,:)) / 2;

        % reflexão (usando valor fixo como no seu código original)
        refl = mid_ + (mid_ - simplexNM(3,:));
        refl(1) = min(max(refl(1), x_min), x_max);
        refl(2) = min(max(refl(2), y_min), y_max);
        fr = nm_val(refl);

        if fr < fvals(2)
            if fr < fvals(1)
                % expansão
                expa = mid_ + 2*(refl - mid_);
                expa(1) = min(max(expa(1), x_min), x_max);
                expa(2) = min(max(expa(2), y_min), y_max);
                fe = nm_val(expa);
                if fe < fr
                    simplexNM(3,:) = expa;
                    fvals(3) = fe;
                else
                    simplexNM(3,:) = refl;
                    fvals(3) = fr;
                end
            else
                simplexNM(3,:) = refl;
                fvals(3) = fr;
            end
        else
            if fr < fvals(3)
                % contração externa
                contr = mid_ + 0.5*(refl - mid_);
            else
                % contração interna
                contr = mid_ + 0.5*(simplexNM(3,:) - mid_);
            end
            contr(1) = min(max(contr(1), x_min), x_max);
            contr(2) = min(max(contr(2), y_min), y_max);
            fc = nm_val(contr);

            if fc < fvals(3)
                simplexNM(3,:) = contr;
                fvals(3) = fc;
            else
                % shrink
                simplexNM(2,:) = simplexNM(1,:) + 0.5*(simplexNM(2,:) - simplexNM(1,:));
                simplexNM(3,:) = simplexNM(1,:) + 0.5*(simplexNM(3,:) - simplexNM(1,:));
                for si = 2:3
                    simplexNM(si,1) = min(max(simplexNM(si,1), x_min), x_max);
                    simplexNM(si,2) = min(max(simplexNM(si,2), y_min), y_max);
                    fvals(si) = nm_val(simplexNM(si,:));
                end
            end
        end

        simplex_history = [simplex_history; simplexNM(1,:)];

        if verbose
            disp(['Iteração ', num2str(iterNM), ' - Melhor valor: ', num2str(-fvals(1))]);
        end
    end

    nm_est = simplexNM(1,:);
end
