
function [x, fval, history] = sbplx(fun, x0, options)
% SBPLX  Subplex algorithm (derivative-free simplex search in subspaces)
%   [x, fval, history] = sbplx(fun, x0, options)

    if nargin < 3, options = struct(); end
    if ~isfield(options, 'maxeval'), options.maxeval = 100; end
    if ~isfield(options, 'tol'), options.tol = 1e-6; end
    if ~isfield(options, 'bounds'), options.bounds = [ -inf(size(x0)); inf(size(x0)) ]; end
    if ~isfield(options, 'prt'), options.prt = 0; end

    n = numel(x0);
    x = x0(:)';  % garantir linha
    lb = options.bounds(1,:);
    ub = options.bounds(2,:);
    maxeval = options.maxeval;
    tol = options.tol;

    h = 0.05 * max(abs(x), 1);  % tamanho inicial do passo
    fx = fun(x);
    neval = 1;
    history = cell(1000, 1);  % pré-alocação
    h_idx = 1;

    while neval < maxeval
        [~, idx] = sort(-abs(h));
        k = 1;

        while k <= n
            m = min(5, n - k + 1);
            idxk = idx(k:k + m - 1);

            simplex = repmat(x, m+1, 1);
            for j = 1:m
                delta = zeros(1,n);
                delta(idxk(j)) = h(idxk(j));
                simplex(j+1,:) = x + delta;
            end
            simplex = min(max(simplex, lb), ub);

            fvals = zeros(m+1,1);
            fvals(1) = fx;
            for j = 2:m+1
                fvals(j) = fun(simplex(j,:));
                neval = neval + 1;
            end

            iter = 0;
            while iter < 20
                [fvals, idxs] = sort(fvals);
                simplex = simplex(idxs,:);
                xbar = mean(simplex(1:m,:), 1);
                xr = 2*xbar - simplex(end,:);
                xr = min(max(xr, lb), ub);
                fr = fun(xr);
                neval = neval + 1;

                if fr < fvals(1)
                    xe = 2*xr - xbar;
                    xe = min(max(xe, lb), ub);
                    fe = fun(xe);
                    neval = neval + 1;
                    if fe < fr
                        simplex(end,:) = xe; fvals(end) = fe;
                    else
                        simplex(end,:) = xr; fvals(end) = fr;
                    end
                elseif fr < fvals(end-1)
                    simplex(end,:) = xr; fvals(end) = fr;
                else
                    xc = 0.5*(simplex(end,:) + xbar);
                    xc = min(max(xc, lb), ub);
                    fc = fun(xc);
                    neval = neval + 1;
                    if fc < fvals(end)
                        simplex(end,:) = xc; fvals(end) = fc;
                    else
                        for j = 2:m+1
                            simplex(j,:) = 0.5*(simplex(j,:) + simplex(1,:));
                            simplex(j,:) = min(max(simplex(j,:), lb), ub);
                            fvals(j) = fun(simplex(j,:));
                            neval = neval + 1;
                        end
                    end
                end

                iter = iter + 1;
                if max(vecnorm(simplex(2:end,:) - simplex(1,:), 2, 2)) < tol
                    break;
                end
            end

            x = simplex(1,:);
            fx = fvals(1);
            history{h_idx} = x;
            h_idx = h_idx + 1;

            k = k + m;
            if neval >= maxeval, break; end
        end

        h = h / 2;
        if max(h) < tol, break; end
    end

    fval = fx;
    history = history(1:h_idx-1);
end
