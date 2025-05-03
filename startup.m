% startup.m – Executado automaticamente ao abrir o projeto
clc;
%clear functions

% ---------------- Parâmetros Gerais ----------------
freq = 15e9;
lambda = (3e8) / freq;
L = 100000;
power = 0.1;
alpha = 2;
P_tx = 0.1;
SNR_dB = 10;
% ---------------- Monte Carlo Simulation  ----------------
MCS = 500;
% --------------- Flags de Plots ---------------------
plt_array = 1;
plt_hiper = 0;
plt_circle = 0;
plt_itersec = 0;
plt_neldermead = 0;
plt_peach = 0;
plt_spectrum = 0;

% --------------- Antenas URA -----------------------
Mx = 4; Mz = 4;
M = Mx * Mz;
d_x = lambda/2;
d_z = lambda/2;
elev = 20;

% --------------- Grade de Busca ---------------------
x_grid = -100:1:100;
y_grid = 10:1:100;
x = [min(x_grid), max(x_grid)];
y = [lambda, max(y_grid)];

x_rand = x_grid(1) + (x_grid(end) - x_grid(1)) * rand;
y_rand = y_grid(1) + (y_grid(end) - y_grid(1)) * rand;
% --------------- Posição do Usuário -----------------

%pos = [x_rand, y_rand, 0];
pos = [30 30 0];
UEs = pos;

% --------------- Parâmetros PEACH -------------------
% n_hiper = 96;
% n_circ = 48;

n_hiper = 12;
n_circ = 6;

% --------------- Parâmetros NM ----------------------
max_iter = 100;
tol = 1e-7;
deltaArea = 1;
numIterNM = 200;

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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Controlar plots logica bool(true=1, false=0)
% 
% plt_array      = 0; % architecture
% plt_hiper      = 0; % PEACH
% plt_circle     = 0; % PEACH
% plt_itersec    = 0; % PEACH
% plt_peach      = 0; % PEACH
% plt_neldermead = 0; % Nelder-Mead
% plt_spectrum   = 1; % pseudospectrum  
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function pseudo_random = u_rand(grid_vector)
%     if min(grid_vector) == 0
%         % Geracao simetrica ao redor de zero
%         pseudo_random = (2 * max(grid_vector)) * rand - max(grid_vector);
%     else
%         % Geracao dentro do intervalo real
%         pseudo_random = min(grid_vector) + (max(grid_vector) - min(grid_vector)) * rand;
%     end
% end