    % startup.m – Executado automaticamente ao abrir o projeto
    clc;
    %clear functions
    
    % ---------------- Parâmetros Gerais ----------------
    freq = 15e9;
    lambda = (3e8) / freq;
    L = 100;
    power = 0.1;
    alpha = 2;
    P_tx = 0.1;
    SNR_dB = 10;
    % ---------------- Monte Carlo Simulation  ----------------
    MCS = 1000;
    % --------------- Flags de Plots ---------------------
    plt_array = 0;
    plt_hiper = 0;
    plt_circle = 0;
    plt_itersec = 0;
    plt_neldermead = 0;
    plt_peach = 0;
    plt_spectrum = 1;
    
    % --------------- Antenas URA -----------------------
    Mx = 8; Mz = 8;
    M = Mx * Mz;
    d_x = lambda/2;
    d_z = lambda/2;
    elev = 20;
    
    % --------------- Grade de Busca (resolução 0.5 m) ---------------------
    % x_grid = -50:1:50;        % reduz região para [−50,50] m e passo de 0.5 m
    % y_grid =  10:1:50;        % idem para [10,50] m
    x_grid = -70:1:70;        % reduz região para [−50,50] m e passo de 0.5 m
    y_grid =  10:1:70;        % idem para [10,50] m
    
    x = [min(x_grid), max(x_grid)];
    y = [min(y_grid), max(y_grid)];
    
    x1 = [-50, 50];
    y1 = [13, 50];
    
    % x_rand = x_grid(1) + (x_grid(end) - x_grid(1)) * rand;
    % y_rand = y_grid(1) + (y_grid(end) - y_grid(1)) * rand;
    
    % x_rand = -50 + (50 - (-50)) * rand;
    % y_rand = 13 + (50 - 10) * rand;
    % --------------- Posição do Usuário -----------------
    
    %pos = [x_rand, y_rand, 0];
    pos = [15 35 0];
    UEs = pos;
    
    % --------------- Parâmetros PEACH -------------------
    n_hiper = 21;
    n_circ = 12;
    
    % --------------- Parâmetros de Refinamento Local Nelder Mead ----------------------
    %max_iter   = 50;            % número de iterações do Hooke–Jeeves/NM
        tol        = 1e-8;          % tolerância para o passo mínimo
        deltaArea  = 0.5;           % passo inicial de 0.5 m, para capturar picos entre nós
        numIterNM  = 15;      % mantém coerência de nomes
    
    % --------------------------- Subplex -------------------------------
    tol_subplex = 1e-5;    % mesma tolerância do Nelder-Mead
    max_eval_subplex = 100; % próximo ao número total de avaliações do Nelder-Mead
    
    
    % 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %                           PARAMETROS
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % freq = 15e9;                      % 15 GHz
    % lambda = (3 * 10^8) / freq;       % comprimento de onda
    % L = 100;                        % numero de amostras temporais
    % power = 0.1;                      % potencia transmitida (W)
    % alpha = 2;                        % expoente path loss
    % P_tx = 0.1;                       % transmission power (w)
    % %N_dBm = -90;                     % noise power in dBm
    % SNR_dB = 19;                       % relacao sinal ruido
    % % ARQUITETURA URA
    % Mx = 8;                           % n antenas na horizontal
    % Mz = 8;                           % n  de antenas na vertical
    % M = Mx * Mz;                      % n  total de antenas
    % d_x = lambda / 2;                 % espacamento antenas eixo x
    % d_z = lambda / 2;                 % espacamento antenas eixo z
    % elev = 20;                        % altura do array
    % % Geracao da grade de busca
    % x_grid = -50:1:50;                % grade de busca em x
    % y_grid = 10:1:50;                 % grade de busca em y
    % % MUSIC
    % max_iter = 10;                    % max interacao do Nelder-Mead
    % tol = 1e-5;                       % tolerancia de erro
    % % PEACH
    % x = [min(x_grid), max(x_grid)];   % limites de busca em x
    % y = [min(y_grid), max(y_grid)];   % limites de busca em y
    % % Posicao da fonte (usuario)
    % x_rand = rand_u(x);               % gerar pos. pseudoaleatoria
    % y_rand = rand_u(y);               % gerar pos. pseudoaleatoria
    % %pos = [rand_u(x) rand_u(y) 0];   % posicao 3D user pseudoaleatoria
    % pos = [30 30 0];                  % posicao determinisca
    % UEs = pos;                        % mesma ref passada ao sinal
    % n_hiper = 816;                     % numero de candidatos
    % n_circ  = 408;                     % numero de candidatos
