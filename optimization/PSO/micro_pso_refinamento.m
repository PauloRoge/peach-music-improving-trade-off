% [est_sbplx, ~] = micro_pso_refinamento(URA, est_peach, ...
%     Un, lambda, ref, x, y, numIterNM);

function [best_pos, history] = micro_pso_refinamento(URA, pos_init, Un, lambda, ref, x_lim, y_lim, max_iter)
    % MICRO_PSO_REFINAMENTO — Refina a estimativa PEACH usando PSO com poucas partículas
    % Entradas:
    %   URA, Un, lambda, ref  -> parâmetros da função MUSIC
    %   pos_init              -> estimativa inicial do PEACH
    %   x_lim, y_lim          -> limites [xmin xmax], [ymin ymax]
    %   max_iter              -> número máximo de iterações
    %
    % Saída:
    %   best_pos              -> estimativa refinada
    %   history               -> histórico de posições

    if nargin < 8
        max_iter = 20;
    end

    % Função objetivo: minimizar o inverso do valor MUSIC
    ps_func = @(xy) -music_eval(xy(1), xy(2), URA, Un, lambda, ref);

    % Parâmetros do PSO reduzido
    n_particles = 3;
    w = 0.6; c1 = 1.5; c2 = 1.5;

    % Inicialização de partículas ao redor da estimativa
    spread = 3; % espalhamento inicial
    particles = repmat(pos_init, n_particles, 1) + spread * (2*rand(n_particles,2) - 1);
    particles(:,1) = min(max(particles(:,1), x_lim(1)), x_lim(2));
    particles(:,2) = min(max(particles(:,2), y_lim(1)), y_lim(2));

    % Velocidade inicial
    velocities = zeros(n_particles, 2);

    % Avaliação inicial
    scores = arrayfun(@(i) ps_func(particles(i,:)), 1:n_particles)';
    pbest = particles;
    pbest_scores = scores;
    [~, gbest_idx] = min(pbest_scores);
    gbest = pbest(gbest_idx,:);

    history = zeros(max_iter, 2);

    for iter = 1:max_iter
        for i = 1:n_particles
            r1 = rand(); r2 = rand();
            velocities(i,:) = w * velocities(i,:) ...
                            + c1 * r1 * (pbest(i,:) - particles(i,:)) ...
                            + c2 * r2 * (gbest - particles(i,:));

            % Atualização da posição
            particles(i,:) = particles(i,:) + velocities(i,:);
            particles(i,1) = min(max(particles(i,1), x_lim(1)), x_lim(2));
            particles(i,2) = min(max(particles(i,2), y_lim(1)), y_lim(2));

            % Avaliação
            score = ps_func(particles(i,:));
            if score < pbest_scores(i)
                pbest(i,:) = particles(i,:);
                pbest_scores(i) = score;
                if score < ps_func(gbest)
                    gbest = particles(i,:);
                end
            end
        end
        history(iter,:) = gbest;
    end
    best_pos = gbest;
end

function val = music_eval(x, y, URA, Un, lambda, ref)
    d_ref = sqrt((ref(1) - x)^2 + (ref(2) - y)^2 + ref(3)^2);
    d_k = sqrt((URA(:,1) - x).^2 + (URA(:,2) - y).^2 + URA(:,3).^2);
    phase_diff = -(2*pi/lambda) * (d_ref - d_k);
    a = exp(1j * phase_diff);
    den = abs(a' * (Un * Un') * a);
    if den < 1e-12
        val = 1e6;
    else
        val = 1 / den;
    end
end
