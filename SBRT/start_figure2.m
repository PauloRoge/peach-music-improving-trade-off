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
SNR_dB = 15;
% ---------------- Monte Carlo Simulation  ----------------
% MCS = 10000;
MCS = 10000;
% --------------- Flags de Plots ---------------------
plt_array = 0;
plt_hiper = 0;
plt_circle = 0;
plt_itersec = 0;
plt_neldermead = 0;
plt_peach = 0;
plt_spectrum = 0;

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
%pos = [10 13 0];
%UEs = pos;

% --------------- Parâmetros PEACH -------------------
n_hiper = 96;
n_circ = 48;

% --------------- Parâmetros de Refinamento Local Nelder Mead ----------------------
%max_iter   = 50;            % número de iterações do Hooke–Jeeves/NM
    tol        = 1e-8;       % tolerância para o passo mínimo
    deltaArea  = 0.5;        % passo inicial de 0.5 m, para capturar picos entre nós
    numIterNM  = 15;         % mantém coerência de nomes

% --------------------------- Subplex -------------------------------
tol_subplex = 1e-5;    % mesma tolerância do Nelder-Mead
max_eval_subplex = 100; % próximo ao número total de avaliações do Nelder-Mead
