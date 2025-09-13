% ---------------- Parâmetros Gerais ----------------
freq   = 15e9;
lambda = (3e8) / freq;   % << CORRIGIDO
L      = 100;
power  = 0.1;
alpha  = 2;
P_tx   = 0.1;

% ---------------- Monte Carlo ----------------------
MCS = 1000;

% ---------------- Flags de Plots -------------------
plt_array      = 0;
plt_hiper      = 0;
plt_circle     = 0;
plt_itersec    = 0;
plt_neldermead = 0;
plt_peach      = 0;
plt_spectrum   = 1;      % << melhor 0 para parfor

% ---------------- Antenas URA ----------------------
Mx = 8; Mz = 8;
M  = Mx * Mz;
d_x = lambda/2;
d_z = lambda/2;
elev = 20;

% --------- Grade e limites -------------------------
x_grid = -70:1:70;
y_grid =  10:1:70;
x = [min(x_grid), max(x_grid)];
y = [min(y_grid), max(y_grid)];

% --------- PEACH e NM ------------------------------
n_hiper   = 96;
n_circ    = 48;
tol       = 1e-8;
deltaArea = 0.5;
numIterNM = 15;
tol_subplex      = 1e-5;
max_eval_subplex = 100;
