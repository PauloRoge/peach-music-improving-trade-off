clc; clear;
SNR_dB =15;
%% ========= Setup (única simulação) ===========================
startup1;   
rng(1234,'twister');% reprodutibilidade da ÚNICA simulação

% Caso o startup_nlos NÃO defina a posição do usuário, 
if ~exist('UEs','var') || isempty(UEs)
    % Gera posição 3D dentro dos limites [x,y], z=0
    UEs = [ ...
        x(1) + (x(2)-x(1))*rand, ...
        y(1) + (y(2)-y(1))*rand, ...
        0];
end

%% ========= Geometria da URA ==================================
[URA, URA_x, URA_z, x_h, x_v, z_h, z_v] = subarrays(Mx, Mz, ...
    d_x, d_z, elev, lambda, plt_array);
URA_y = zeros(size(URA_x));
ref   = URA(1,:); % elemento de referência

%% ========= Geração do sinal (NLoS multipercurso) =============
K_dB = 3;  % Rician K (em dB) – aqui 0 dB: LoS e NLoS

scatterer_pos = [  % (Lx3) posições dos espalhadores
    -35, 15, 2;
     -8, 10, 1.2;
     20,  5, 2.5
];

Gamma = [ ...      % coeficientes de reflexão (módulo e fase)
    0.5*exp(1j*pi/3);
    0.3*exp(1j*1.1*pi);
    0.4*exp(1j*0.2*pi)
];

% Saída: Yh, Yv (subarranjos) e Y (M x L)
[Yh, Yv, Y] = signals_nlos_multi(UEs, URA, lambda, L, ...
    alpha, SNR_dB, P_tx, Mx, Mz, K_dB, scatterer_pos, Gamma);

%% ========= Pseudoespectro (para inspeção/depuração opcional) =
% Função de steering do array (3D). 
responsearray = @(xq, yq, zq) steering_vector(Mx, Mz, ...
    elev, d_x, d_z, lambda, xq, yq, zq);
[Pmusic] = pseudospectrum(responsearray, Y, L);

%% ========= Estimativa PEACH-Golden (única execução) ==========
% Atenção à ordem dos argumentos: use x,y (limites) 

[Un_h, Un_v, est_peach] = peach_golden( ...
    Yh, Yv, L, x, n_hiper, x_h, z_h, x_v, z_v, ref, lambda, y, n_circ, UEs);

fprintf('\nEstimativa PEACH-Golden: (x=%.2f, y=%.2f)\n', ...
    est_peach(1), est_peach(2));

%% ========= Subespaço de ruído (MUSIC completo para NM) =======
Cov = (Y*Y')/L;
[V,D] = eig(Cov);
[~, ord] = sort(diag(D), 'descend');
V = V(:, ord);

estimated_sources = 1;                  % cenário 1 usuário
Un = V(:, estimated_sources+1:end);     % subespaço de ruído

%% ========= Refinamento local (Nelder–Mead) ===================
% deltaArea, numIterNM, tolerâncias vêm do startup_nlos
[nm_est, history] = nelder_mead(URA, est_peach, Un, ...
    lambda, ref, deltaArea, numIterNM, 1e-6, true, x, y); 

%% ========= Métricas (erros e CRB) ============================
erro_peach  = norm(est_peach(1:2) - UEs(1:2));
erro_nelder = norm(nm_est(1:2)     - UEs(1:2));

fprintf('Posição real do usuário: (x=%.2f, y=%.2f)\n', ...
    UEs(1), UEs(2));
fprintf('Erro PEACH-Golden: %.4f m\n', erro_peach);
fprintf('Erro Nelder–Mead:  %.4f m\n', erro_nelder);

crb_eucl = crb_vetorizado(L, URA, UEs, lambda, P_tx, SNR_dB);
fprintf('CRB (Euclidiano):  %.4f m\n', crb_eucl);
printResults;
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       FUNÇÕES LOCAIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = steering_vector(Mx, Mz, elev, d_x, d_z, lambda ...
    , x, y, z)
    x_pos = (0:Mx-1) * d_x;
    z_pos = (0:Mz-1) * d_z + elev;

    a = zeros(Mx * Mz, 1);
    idx = 1;

    % Referência (1,1)
    x_ref = x_pos(1);
    z_ref = z_pos(1);
    d_ref = sqrt((x - x_ref)^2 + y^2 + (z - z_ref)^2);

    for i = 1:Mx
        for j = 1:Mz
            x_k = x_pos(i);
            z_k = z_pos(j);
            d_k = sqrt((x - x_k)^2 + y^2 + (z - z_k)^2);
            phase_k = -(2*pi/lambda) * (d_ref - d_k);
            a(idx) = exp(1j * phase_k);
            idx = idx + 1;
        end
    end
end

function [Pmusic, Un] = pseudospectrum(responsearray, Y, ...
    snapshots)
    R = (Y * Y') / snapshots;
    [V,D] = eig(R);
    [~, iord] = sort(diag(D), 'descend');
    V = V(:, iord);

    estimated_sources = 1;
    Un = V(:, estimated_sources+1:end);

    Pmusic = @(pos) pmusic(responsearray, pos, Un);
end

function val = pmusic(responsearray, pos, Un)
    % pos = [x y z]
    a = responsearray(pos(1), pos(2), pos(3));
    % Projeção no subespaço de ruído
    ProjN = Un * (Un');
    denom = abs(a' * ProjN * a);
    % Proteção numérica mínima
    if denom < eps, denom = eps; end
    val = 1 / denom;
end
