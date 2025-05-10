% parâmetros figura (1)
clc;
% ---------------- Parâmetros Gerais ----------------
freq = 15e9;
lambda = (3e8) / freq;
L = 100;
power = 1;
alpha = 2;
P_tx = 0.1;
SNR_dB = 15;
% ---------------- Monte Carlo Simulation  ----------------
MCS = 2000; 
% --------------- Antenas URA -----------------------
Mx = 4; Mz = 4;
M = Mx * Mz;
d_x = lambda/2;
d_z = lambda/2;
elev = 15;
% --------------- Grade de Busca ---------------------
x_grid = -70:1:70;
y_grid = 10:1:70;
x = [min(x_grid), max(x_grid)];
y = [lambda, max(y_grid)];

%x_rand = x_grid(1) + (x_grid(end) - x_grid(1)) * rand;
%y_rand = y_grid(1) + (y_grid(end) - y_grid(1)) * rand;
%pos = [x_rand, y_rand, 0];

% --------------- Posição do Usuário -----------------
pos = [10 30 0];
UEs = pos;

% --------------- Parâmetros PEACH -------------------
n_hiper = 200;
n_circ = 100;

% --------------- Parâmetros NM ----------------------
max_iter = 100;
tol = 1e-6;
deltaArea = 3;
numIterNM = 100;



% % parâmetros figura (1)
% clc;
% % ---------------- Parâmetros Gerais ----------------
% freq = 15e9;
% lambda = (3e8) / freq;
% L = 100000;
% power = 0.1;
% alpha = 2;
% P_tx = 0.1;
% SNR_dB = 10;
% % ---------------- Monte Carlo Simulation  ----------------
% MCS = 10; 
% % --------------- Antenas URA -----------------------
% Mx = 4; Mz = 4;
% M = Mx * Mz;
% d_x = lambda/2;
% d_z = lambda/2;
% elev = 20;
% % --------------- Grade de Busca ---------------------
% x_grid = -100:1:100;
% y_grid = 10:1:100;
% x = [min(x_grid), max(x_grid)];
% y = [lambda, max(y_grid)];
% 
% %x_rand = x_grid(1) + (x_grid(end) - x_grid(1)) * rand;
% %y_rand = y_grid(1) + (y_grid(end) - y_grid(1)) * rand;
% %pos = [x_rand, y_rand, 0];
% 
% % --------------- Posição do Usuário -----------------
% pos = [10 30 0];
% UEs = pos;
% 
% % --------------- Parâmetros PEACH -------------------
% n_hiper = 96;
% n_circ = 48;
% 
% % --------------- Parâmetros NM ----------------------
% max_iter = 100;
% tol = 1e-6;
% deltaArea = 3;
% numIterNM = 100;