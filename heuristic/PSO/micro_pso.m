% % -------- Micro-PSO em vez do PEACH -----------------
% nPart = 6;        % partículas
% nIter = 20;       % iterações
% Pfun  = @(x,y) Pmusic([x, y, 0]);   % wrapper 2D
% 
% [pso_est, pswarm] = micro_pso(Pfun, x, y, nPart, nIter);
% fprintf('\nPosição Micro-PSO: (%.2f, %.2f)', pso_est(1), pso_est(2));

function [pso_est, swarm_hist] = micro_pso(Pmusic, x_lim, y_lim, nPart, nIter)
% MICRO_PSO  PSO compacto (2 D) para estimar [x,y] do pico do pseudo-espectro
%
%   Entrada
%     Pmusic  : handle @(x,y)→valor  (quanto maior, melhor)
%     x_lim   : [xmin xmax]
%     y_lim   : [ymin ymax]
%     nPart   : nº de partículas   (4–8 é suficiente)
%     nIter   : nº de iterações    (15–30 tipicamente)
%
%   Saída
%     pso_est    : melhor posição encontrada  [x y]
%     swarm_hist : trajetória das partículas (cell{iter}→nPart×2)

% parâmetros clássicos reduzidos
w  = 0.6;      % inércia
c1 = 1.7;      % atração ao pbest
c2 = 1.7;      % atração ao gbest

% inicialização aleatória uniforme
X = [ x_lim(1) + (x_lim(2)-x_lim(1))*rand(nPart,1), ...
      y_lim(1) + (y_lim(2)-y_lim(1))*rand(nPart,1) ];
V = zeros(nPart,2);

Pbest   = X;
PbestVal= arrayfun(@(i) Pmusic(X(i,1),X(i,2)), 1:nPart)';
[gbestVal, idx] = max(PbestVal);   gbest = Pbest(idx,:);

swarm_hist = cell(nIter,1);

for it = 1:nIter
    r1 = rand(nPart,2); r2 = rand(nPart,2);
    V  = w*V + c1*r1.*(Pbest - X) + c2*r2.*(gbest - X);
    X  = X + V;
    % limites
    X(:,1) = min(max(X(:,1), x_lim(1)), x_lim(2));
    X(:,2) = min(max(X(:,2), y_lim(1)), y_lim(2));

    % avaliação
    f = arrayfun(@(i) Pmusic(X(i,1),X(i,2)), 1:nPart)';
    improve            = f > PbestVal;
    Pbest(improve,:)   = X(improve,:);
    PbestVal(improve)  = f(improve);
    [gbestVal, idx]    = max(PbestVal);
    gbest              = Pbest(idx,:);

    swarm_hist{it} = X;   % para plot/depuração
end

pso_est = gbest;
end
