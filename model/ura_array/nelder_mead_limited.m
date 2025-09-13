function [nm_est, simplex_history] = nelder_mead_limited(URA, pos_peach, Un, lambda, ref, delta, max_iter, tol, verbose, x, y)

    if nargin < 13, tol = 1e-6; end
    if nargin < 14, verbose = false; end

    % Função objetivo
    ps_func = @(x,y) music(x, y, URA, Un, lambda, ref);
    nm_val  = @(xy) -ps_func(xy(1), xy(2));

    % Inicializa simplex
    simplexNM = [
        pos_peach(1), pos_peach(2);
        pos_peach(1) + delta, pos_peach(2);
        pos_peach(1), pos_peach(2) + delta
    ];

    % Garante que o simplex inicial fique nos limites
    simplexNM(:,1) = min(max(simplexNM(:,1), x(1)), x(2));
    simplexNM(:,2) = min(max(simplexNM(:,2), y(1)), y(2));

    fvals = zeros(3,1);
    for si = 1:3
        fvals(si) = nm_val(simplexNM(si,:));
    end

    simplex_history = cell(max_iter+1,1);
    simplex_history{1} = simplexNM;

    for iterNM = 1:max_iter
        [fvals, idxSort] = sort(fvals);
        simplexNM = simplexNM(idxSort,:);
        mid_ = (simplexNM(1,:) + simplexNM(2,:)) / 2;

        refl = 2*mid_ - simplexNM(3,:);

        %limite
        refl(1) = min(max(refl(1), x(1)), x(2));
        refl(2) = min(max(refl(2), y(1)), y(2));
        
        fr = nm_val(refl);

        if fr < fvals(1)
            expa = 2*refl - mid_;

            %limite
            expa(1) = min(max(expa(1), x(1)), x(2));
            expa(2) = min(max(expa(2), y(1)), y(2));
            
            fe = nm_val(expa);
            if fe < fr
                simplexNM(3,:) = expa; fvals(3) = fe;
            else
                simplexNM(3,:) = refl; fvals(3) = fr;
            end
        elseif fr < fvals(2)
            simplexNM(3,:) = refl; fvals(3) = fr;
        else
            contr = (simplexNM(3,:) + mid_) / 2;

            %limite
            contr(1) = min(max(contr(1), x(1)), x(2));
            contr(2) = min(max(contr(2), y(1)), y(2));
            
            fc = nm_val(contr);
            if fc < fvals(3)
                simplexNM(3,:) = contr; fvals(3) = fc;
            else
                % Shrink (contração)
                simplexNM(2,:) = (simplexNM(1,:) + simplexNM(2,:)) / 2;
                simplexNM(3,:) = (simplexNM(1,:) + simplexNM(3,:)) / 2;
                
                % limite
                simplexNM(2,1) = min(max(simplexNM(2,1), x(1)), x(2));
                simplexNM(2,2) = min(max(simplexNM(2,2), y(1)), y(2));
                simplexNM(3,1) = min(max(simplexNM(3,1), x(1)), x(2));
                simplexNM(3,2) = min(max(simplexNM(3,2), y(1)), y(2));

                fvals(2) = nm_val(simplexNM(2,:));
                fvals(3) = nm_val(simplexNM(3,:));
            end
        end

        simplex_history{iterNM+1} = simplexNM;

        if norm(simplexNM(3,:) - simplexNM(1,:)) < tol
            if verbose
                fprintf('Convergência atingida na iteração %d\n', iterNM);
            end
            break;
        end
    end

    nm_est = simplexNM(1,:);
    simplex_history = simplex_history(~cellfun(@isempty, simplex_history));
end

function valPS = music(x, y, URA, Un, lambda, ref)
    % music: Retorna o valor do pseudo-espectro MUSIC em (x,y)
    % URA: matriz M×3 com as posições [x_k, y_k, z_k]

    d_ref = sqrt((ref(1) - x)^2 + (ref(2) - y)^2 + ref(3)^2);

    % Vetor de distâncias d_km
    d_km = sqrt((URA(:,1) - x).^2 + (URA(:,2) - y).^2 + URA(:,3).^2);

    % Vetor de steering (vetorizado)
    phase_diff = -(2*pi/lambda) * (d_ref - d_km);
    a = exp(1j * phase_diff);

    den = abs(a' * (Un * Un') * a);
    if den < 1e-12
        valPS = 1e6;
    else
        valPS = 1 / den;
    end
end
