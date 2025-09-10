% [nm_est, simplex_history] = nelder_mead_test(URA, est_peach, Un, lambda, ref, ...
%                                               deltaArea, numIterNM, tol, 200, true, x, y);

function [nm_est, simplex_history, iter] = nelder_mead_test(URA, pos_peach, Un, lambda, ref, ...
                                                 delta, max_iter, tol, budget,verbose, x, y)
% NELDER_MEAD com contagem de avaliações e print do critério de parada.

    ps_fun = @(xy) -music(xy(1), xy(2), URA, Un, lambda, ref);
    neval = 0;

    % Inicialização do simplex
    simplex = [
        pos_peach;
        pos_peach + [delta, 0];
        pos_peach + [0, delta]
    ];
    simplex(:,1) = min(max(simplex(:,1), x(1)), x(2));
    simplex(:,2) = min(max(simplex(:,2), y(1)), y(2));

    % Avaliação inicial
    fvals = zeros(3,1);
    for i = 1:3
        fvals(i) = ps_fun(simplex(i,:));
        neval = neval + 1;
        if neval >= budget
            %fprintf('nm parou por atingir budget de %d avaliações (inicialização)\n', neval);
            nm_est = simplex(i,:);
            simplex_history = simplex(1:i,:);
            return;
        end
    end

    simplex_history = cell(max_iter+1,1);
    simplex_history{1} = simplex;

    for iter = 1:max_iter
        [fvals, idx] = sort(fvals);
        simplex = simplex(idx,:);
        xbar = mean(simplex(1:2,:),1);

        % Reflexão
        xr = 2*xbar - simplex(3,:);
        xr = min(max(xr, x(1)), x(2));
        xr(2) = min(max(xr(2), y(1)), y(2));
        fr = ps_fun(xr); neval = neval + 1;
        if neval >= budget
            %fprintf('nm parou por atingir budget de %d avaliações (reflexão)\n', neval);
            nm_est = simplex(1,:);
            simplex_history(iter+1) = {simplex(1,:)};
            return;
        end

        if fr < fvals(1)
            xe = 2*xr - xbar;
            xe(1) = min(max(xe(1), x(1)), x(2));
            xe(2) = min(max(xe(2), y(1)), y(2));
            fe = ps_fun(xe); neval = neval + 1;
            if fe < fr
                simplex(3,:) = xe; fvals(3) = fe;
            else
                simplex(3,:) = xr; fvals(3) = fr;
            end
        elseif fr < fvals(2)
            simplex(3,:) = xr; fvals(3) = fr;
        else
            xc = 0.5*(simplex(3,:) + xbar);
            xc(1) = min(max(xc(1), x(1)), x(2));
            xc(2) = min(max(xc(2), y(1)), y(2));
            fc = ps_fun(xc); neval = neval + 1;
            if neval >= budget
                %fprintf('nm parou por atingir budget de %d avaliações (contração)\n', neval);
                nm_est = simplex(1,:);
                simplex_history(iter+1) = {simplex(1,:)};
                return;
            end
            if fc < fvals(3)
                simplex(3,:) = xc; fvals(3) = fc;
            else
                for j = 2:3
                    simplex(j,:) = simplex(1,:) + 0.5*(simplex(j,:) - simplex(1,:));
                    simplex(j,1) = min(max(simplex(j,1), x(1)), x(2));
                    simplex(j,2) = min(max(simplex(j,2), y(1)), y(2));
                    fvals(j) = ps_fun(simplex(j,:));
                    neval = neval + 1;
                    if neval >= budget
                        %fprintf('nm parou por atingir budget de %d avaliações (shrink)\n', neval);
                        nm_est = simplex(1,:);
                        simplex_history(iter+1) = {simplex(1,:)};
                        return;
                    end
                end
            end
        end

        simplex_history{iter+1} = simplex(1,:);
        if max(vecnorm(simplex(2:3,:) - simplex(1,:), 2, 2)) < tol
            %fprintf('nm convergiu por tolerância após %d iterações e %d avaliações\n', iter, neval);
            break;
        end
    end

    % Se finalizou por max_iter
    if iter == max_iter
        %fprintf('nm parou por atingir o número máximo de iterações (%d).\n', max_iter);
    end

    nm_est = simplex(1,:);
    simplex_history = simplex_history(~cellfun('isempty',simplex_history));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCAO LOCAL - MUSIC DINAMIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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