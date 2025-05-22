function [x, fval, history, neval] = sbplx_test(fun, x0, options)
% SBPLX  Subplex algorithm com prints informando critério de parada.
%   [x, fval, history, neval] = sbplx(fun, x0, options)

    if nargin < 3, options = struct(); end
    if ~isfield(options, 'maxeval'), options.maxeval = 100; end
    if ~isfield(options, 'tol'),    options.tol     = 1e-6; end
    if ~isfield(options, 'bounds'), options.bounds  = [ -inf(size(x0)); inf(size(x0)) ]; end
    if ~isfield(options, 'prt'),    options.prt     = 0; end

    n = numel(x0);
    x = x0(:)';
    lb = options.bounds(1,:);
    ub = options.bounds(2,:);
    maxeval = options.maxeval;
    tol     = options.tol;

    h = 0.05 * max(abs(x), 1);
    fx = fun(x);
    neval   = 1;
    history = cell(1000,1);
    h_idx   = 1;

    % Loop principal
    while neval < maxeval
        [~, idx] = sort(-abs(h));
        k = 1;
        while k <= n && neval < maxeval
            m = min(5, n - k + 1);
            idxk = idx(k:k+m-1);

            % Monta simplex
            simplex = repmat(x, m+1, 1);
            for j = 1:m
                delta = zeros(1,n);
                delta(idxk(j)) = h(idxk(j));
                simplex(j+1,:) = x + delta;
            end
            simplex = min(max(simplex, lb), ub);

            % Avalia vértices
            fvals = zeros(m+1,1);
            fvals(1) = fx;
            for j = 2:m+1
                fvals(j) = fun(simplex(j,:));
                neval = neval + 1;
                if neval >= maxeval, break; end
            end

            % Sub-simplex NM
            iter = 0;
            while iter < 20 && neval < maxeval
                [fvals, ids] = sort(fvals);
                simplex = simplex(ids,:);
                xbar = mean(simplex(1:m,:),1);

                % Reflexão
                xr = 2*xbar - simplex(end,:);
                xr = min(max(xr, lb), ub);
                fr = fun(xr); neval = neval + 1;
                if fr < fvals(1)
                    % Expansão
                    xe = 2*xr - xbar;
                    xe = min(max(xe, lb), ub);
                    fe = fun(xe); neval = neval + 1;
                    if fe < fr
                        simplex(end,:) = xe; fvals(end) = fe;
                    else
                        simplex(end,:) = xr; fvals(end) = fr;
                    end
                elseif fr < fvals(end-1)
                    simplex(end,:) = xr; fvals(end) = fr;
                else
                    % Contração
                    xc = 0.5*(simplex(end,:) + xbar);
                    xc = min(max(xc, lb), ub);
                    fc = fun(xc); neval = neval + 1;
                    if fc < fvals(end)
                        simplex(end,:) = xc; fvals(end) = fc;
                    else
                        % Shrink
                        for j = 2:m+1
                            simplex(j,:) = 0.5*(simplex(j,:) + simplex(1,:));
                            simplex(j,:) = min(max(simplex(j,:), lb), ub);
                            fvals(j) = fun(simplex(j,:));
                            neval = neval + 1;
                            if neval >= maxeval, break; end
                        end
                    end
                end

                iter = iter + 1;
                % Convergência por tolerância no simplex
                if max(vecnorm(simplex(2:end,:) - simplex(1,:),2,2)) < tol
                    %fprintf('\nsbplx convergiu por tolerância no simplex após %d avaliações\n', neval);
                    break;
                end
            end

            % Atualiza ponto corrente
            x  = simplex(1,:);
            fx = fvals(1);
            history{h_idx} = x;
            h_idx = h_idx + 1;
            k = k + m;
        end

        % Reduz passo e verifica tolerância global
        h = h/2;
        if max(h) < tol
            %fprintf('\nsbplx convergiu por tamanho de passo < tol após %d avaliações\n', neval);
            break;
        end
    end

    % Se esgotou maxeval
    if neval >= maxeval
        %fprintf('\nsbplx terminou por atingir maxeval (%d avaliações)\n', neval);
    end

    fval   = fx;
    history = history(1:h_idx-1);
end
