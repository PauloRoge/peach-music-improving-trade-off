%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PEACH -> Nelder-Mead (SEM busca coarse), c/ pseudo-espectro MUSIC
clear; close all; clc;
tic;

%% Parâmetros Principais
SNR_dB     = 20;
z_antenna  = 0;  % altura do array
L          = 10;  % snapshots
user_pos   = [-25/sqrt(2), 42/sqrt(2), 0]; % Exemplo
dd         = 30;  % distância nominal
f_c        = 15e9;
lambda     = 3e8/f_c;
d          = lambda/2;
Nx         = 8;   % arranjo Nx x Ny
Ny         = 8;
N          = Nx*Ny;
N_hiper    = 1000; 
N_circ     = 500;
x_lim      = 1.5*dd;
y_lim      = 1.5*dd;

%% Parâmetros do refinamento Nelder-Mead
deltaArea  = 5;   % "Raio" para inicializar o simplex
numIterNM  = 100; % Nº máximo de iterações Nelder-Mead

%% 1) Geometria do URA c/ z=20m
x_off = (Nx-1)*d/2;
z_off = (Ny-1)*d/2;
elements_x_a = zeros(Nx,Ny);
elements_y_a = zeros(Nx,Ny);
elements_z_a = zeros(Nx,Ny);
for ix_ = 1:Nx
    for iy_ = 1:Ny
        elements_x_a(ix_, iy_) = (ix_-1)*d - x_off;
        elements_y_a(ix_, iy_) = 0;
        elements_z_a(ix_, iy_) = (iy_-1)*d - z_off + z_antenna;
    end
end
x_ref = elements_x_a(1,1);
y_ref = elements_y_a(1,1);
z_ref = elements_z_a(1,1);

%% 2) Geração do Sinal (1 realização)
SNR_lin = 10^(SNR_dB/10);
X_2D = zeros(N,L);
X_h  = zeros(Nx,L);  % subarray horizontal (iy=1)
X_v  = zeros(Ny,L);  % subarray vertical   (ix=1)

total_rx_power = 0; % acumula potência recebida (sinal puro, sem ruído)

for t=1:L
    idx_temp = 1;
    d_k1 = sqrt( (x_ref-user_pos(1))^2 + ...
                 (y_ref-user_pos(2))^2 + ...
                 (z_ref-user_pos(3))^2 );

    for ix2=1:Nx
        for iy2=1:Ny
            x_km = elements_x_a(ix2, iy2);
            y_km = elements_y_a(ix2, iy2);
            z_km = elements_z_a(ix2, iy2);

            d_km = sqrt( (x_km - user_pos(1))^2 + ...
                         (y_km - user_pos(2))^2 + ...
                         (z_km - user_pos(3))^2 );

            beta_ = (lambda/(4*pi*d_km))^2;
            phase_ = -(2*pi/lambda)*(d_k1 - d_km);
            val_ = sqrt(beta_)*exp(1j*phase_);

            X_2D(idx_temp,t) = X_2D(idx_temp,t) + val_;
            if iy2==1, X_h(ix2,t) = X_h(ix2,t) + val_; end
            if ix2==1, X_v(iy2,t) = X_v(iy2,t) + val_; end
            idx_temp = idx_temp+1;

            if t==1
                total_rx_power = total_rx_power + beta_;
            end
        end
    end
end

% Potência do ruído necessária (em função da potência do sinal recebido)
noise_power = total_rx_power / SNR_lin;

% Adiciona ruído ao sinal
for t=1:L
    X_2D(:,t) = X_2D(:,t) + sqrt(noise_power/2)*(randn(N,1)+1j*randn(N,1));
    X_h(:,t)  = X_h(:,t)  + sqrt(noise_power/2)*(randn(Nx,1)+1j*randn(Nx,1));
    X_v(:,t)  = X_v(:,t)  + sqrt(noise_power/2)*(randn(Ny,1)+1j*randn(Ny,1));
end

% SNR efetiva realista (definida como sinal puro vs ruído gerado)
SNR_efetiva_dB = 10*log10(total_rx_power / noise_power);

fprintf('SNR desejada (dB): %.2f\n', SNR_dB);
fprintf('SNR efetiva simulada (dB): %.2f\n', SNR_efetiva_dB);


%% 3) Checando SNR efetiva da simulação

% Potência média do sinal recebido (sem ruído)
signal_power = mean(abs(X_2D(:)).^2);

% Potência do ruído que você adicionou
noise_power = total_rx_power / SNR_lin;

% SNR efetiva da simulação
snr_efetiva_dB = 10*log10(signal_power / noise_power);

fprintf('SNR desejada (dB): %.2f\n', SNR_dB);
fprintf('SNR efetiva simulada (dB): %.2f\n', snr_efetiva_dB);
    

%% 3) Covariância & Subespaços
Rxx_2D = (1/L)*(X_2D*X_2D');
[U2D,~,~] = svd(Rxx_2D);
Un_2D = U2D(:,2:end);

Rxx_h = (1/L)*(X_h*X_h');
[Uh,~,~] = svd(Rxx_h);  Un_h = Uh(:,2:end);

Rxx_v = (1/L)*(X_v*X_v');
[Uv,~,~] = svd(Rxx_v);  Un_v = Uv(:,2:end);

%% 4) PEACH 
% (a) Hiperbola (eixo x)
x_candidates = linspace(-x_lim, x_lim, N_hiper);
Pmusic_hvals = zeros(size(x_candidates));
for ixx=1:length(x_candidates)
    xx_cand = x_candidates(ixx);
    a_h_ = zeros(Nx,1);
    d_ref_ = sqrt( (x_ref - xx_cand)^2 + (z_ref - 0)^2 );
    for ix2=1:Nx
        x_km = elements_x_a(ix2,1);
        z_km = elements_z_a(ix2,1);
        d_km_ = sqrt((x_km - xx_cand)^2 + (z_km - 0)^2);
        ph = -(2*pi/lambda)*(d_ref_ - d_km_);
        a_h_(ix2) = exp(1j*ph);
    end
    Pmusic_hvals(ixx) = 1/abs( a_h_'*(Un_h*Un_h')*a_h_ );
end
[~, idx_mh] = max(Pmusic_hvals);
x_peak_h = x_candidates(idx_mh);

% Delta_est
F1x = elements_x_a(1,1);
F1z = elements_z_a(1,1);
F2x = elements_x_a(Nx,1);
F2z = elements_z_a(Nx,1);
dF1 = sqrt( (x_peak_h - F1x)^2 + F1z^2 );
dF2 = sqrt( (x_peak_h - F2x)^2 + F2z^2 );
Delta_est = dF1 - dF2;

% (b) Círculo (eixo y)
y_candidates = linspace(0, y_lim, N_circ);
Pmusic_vvals = zeros(size(y_candidates));
for iyy=1:length(y_candidates)
    yy_cand = y_candidates(iyy);
    a_v_ = zeros(Ny,1);
    d_ref_ = sqrt((x_ref-0)^2 + (y_ref - yy_cand)^2 + (z_ref - 0)^2);
    for iy2=1:Ny
        x_km = elements_x_a(1, iy2);
        z_km = elements_z_a(1, iy2);
        d_km_ = sqrt((x_km-0)^2 + (yy_cand-0)^2 + z_km^2);
        ph2 = -(2*pi/lambda)*( d_ref_ - d_km_ );
        a_v_(iy2) = exp(1j*ph2);
    end
    Pmusic_vvals(iyy) = 1/abs( a_v_'*(Un_v*Un_v')*a_v_ );
end
[~, idx_mv] = max(Pmusic_vvals);
y_peak_v = y_candidates(idx_mv);
R_est = y_peak_v;

% Resolver Interseção => (xPEACH, yPEACH)
syms xs ys real
expr_d1 = sqrt( (xs - F1x)^2 + ys^2 + F1z^2 );
expr_d2 = sqrt( (xs - F2x)^2 + ys^2 + F2z^2 );
eq_hip  = expr_d1 - expr_d2 - Delta_est;
eq_circ = xs^2 + ys^2 - R_est^2;
sol_    = solve([eq_hip, eq_circ],[xs, ys],'Real',true);
x_sol_  = double(sol_.xs);
y_sol_  = double(sol_.ys);

if isempty(x_sol_)
    xPEACH = NaN;  yPEACH=NaN;
else
    distBest = Inf;
    for kk=1:length(x_sol_)
        candPos = [ x_sol_(kk), y_sol_(kk) ];
        dd_ = norm( candPos - user_pos(1:2) );
        if dd_ < distBest
            distBest = dd_;
            xPEACH = candPos(1);
            yPEACH = candPos(2);
        end
    end
end
fprintf('PEACH final = (%.3f, %.3f)\n', xPEACH, yPEACH);

%% 5) Nelder-Mead (refinando c/ Nelder-Mead -> sem coarse local)
% A ideia é maximizar pseudo-espectro => definimos ps_func e nm_val
ps_func = @(xx,yy) computePseudo_2D(xx, yy, ...
    elements_x_a, elements_y_a, elements_z_a, ...
    Un_2D, lambda, x_ref, y_ref, z_ref);
nm_val  = @(xy) -ps_func(xy(1), xy(2));  % - para "maximizar"

% Inicia simplex a partir do PEACH
x_min_loc = xPEACH - deltaArea;
y_min_loc = yPEACH - deltaArea;
simplexNM = [
  xPEACH,        yPEACH;
  xPEACH+1,      yPEACH;     % Deslocamento ex.: 1 m
  xPEACH,        yPEACH+1
];
fvals = zeros(3,1);
for si=1:3
    fvals(si) = nm_val(simplexNM(si,:));
end

for iterNM=1:numIterNM
    % Ordena
    [fvals, idxSort] = sort(fvals);
    simplexNM = simplexNM(idxSort,:);
    mid_ = (simplexNM(1,:) + simplexNM(2,:))/2;
    refl = 2*mid_ - simplexNM(3,:);
    fr   = nm_val(refl);

    if fr < fvals(1)
        % expansion
        expa = 2*refl - mid_;
        fe   = nm_val(expa);
        if fe < fr
            simplexNM(3,:) = expa;  fvals(3)=fe;
        else
            simplexNM(3,:) = refl;  fvals(3)=fr;
        end
    elseif fr < fvals(2)
        simplexNM(3,:) = refl;  fvals(3)=fr;
    else
        contr = (simplexNM(3,:) + mid_)/2;
        fc    = nm_val(contr);
        if fc < fvals(3)
            simplexNM(3,:) = contr; fvals(3)=fc;
        else
            % shrink
            simplexNM(2,:) = (simplexNM(1,:) + simplexNM(2,:))/2;
            simplexNM(3,:) = (simplexNM(1,:) + simplexNM(3,:))/2;
            fvals(2) = nm_val(simplexNM(2,:));
            fvals(3) = nm_val(simplexNM(3,:));
        end
    end

    if norm(simplexNM(3,:) - simplexNM(1,:))<1e-6
        break;
    end
end
nm_est = simplexNM(1,:);
fprintf('Nelder-Mead final = (%.3f, %.3f)\n', nm_est(1), nm_est(2));
fprintf('User real = (%.3f, %.3f)\n', user_pos(1), user_pos(2));

%% 6) Mapa de Calor: pseudo-espectro no range x=-50..50, y=0..50
xx_grid = linspace(-50, 50, 100);
yy_grid = linspace(0, 50, 100);
HeatMap = zeros(length(xx_grid), length(yy_grid));
for ixg=1:length(xx_grid)
    for iyg=1:length(yy_grid)
        HeatMap(ixg, iyg) = ps_func(xx_grid(ixg), yy_grid(iyg));
    end
end

HeatMap = zeros(length(yy_grid), length(xx_grid));
for iy=1:length(yy_grid)
    for ix=1:length(xx_grid)
        HeatMap(iy, ix) = ps_func(xx_grid(ix), yy_grid(iy));
    end
end

%% 7) Plot Mapa (Heatmap) + Localização
%arrume esse plot
figure('Name','Pseudo-espectro (Full 8x8) - sem busca coarse','Color',[1 1 1]);
imagesc(yy_grid, xx_grid, HeatMap); 
axis xy; hold on; colorbar;
xlabel('x (m)');
ylabel('y (m)');
title(sprintf('Pseudo-Espectro (8x8) - SNR=%d dB, z=%d m', SNR_dB, z_antenna));

% Reordene os pares (x, y) no plot:
plot(user_pos(2), user_pos(1), 'wp','MarkerSize',10,...
     'MarkerFaceColor','k','DisplayName','User Real Location');
plot(yPEACH, xPEACH, 'wo','MarkerSize',8,...
     'MarkerFaceColor','g','DisplayName','PEACH Estimation');
plot(nm_est(2), nm_est(1), 'ws','MarkerSize',8,...
     'MarkerFaceColor','r','DisplayName','PEACH-Nelde Estimation');

legend('Location','best');

%% 8) CRB => valor final p/ exibir
[Jxx_c, Jyy_c, Jxy_c] = precomputeCRBterms(elements_x_a, elements_y_a,...
    elements_z_a, user_pos, lambda);
noise_power = total_rx_power / SNR_lin;
P_norm = L*(lambda/(sqrt(2)*pi*sqrt(noise_power)))^2;
Jxx_s = P_norm*Jxx_c; 
Jyy_s = P_norm*Jyy_c;
Jxy_s = P_norm*Jxy_c;
den_x = (Jxx_s - (Jxy_s^2/Jyy_s));
den_y = (Jyy_s - (Jxy_s^2/Jxx_s));
crb_x = 1/den_x;  crb_y = 1/den_y;
crb_eucl = sqrt(crb_x + crb_y);
fprintf('\nCRB (Eucl) = %.4f m\n', crb_eucl);

%% Cálculo e impressão dos erros euclidianos (PEACH e Nelder-Mead)
erro_peach = norm([xPEACH, yPEACH] - user_pos(1:2));
erro_nelder = norm(nm_est - user_pos(1:2));
fprintf('Erro Eucl. PEACH = %.3f m | Erro Eucl. Nelder-Mead = %.3f m\n', ...
    erro_peach, erro_nelder);


elapsedTime = toc;
fprintf('\nConcluído. Tempo total = %.2f s.\n', elapsedTime);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Funções Auxiliares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function valPS = computePseudo_2D(xx, yy, elements_x, elements_y, ...
    elements_z, Un_2D, lambda, x_ref, y_ref, z_ref)
% computePseudo_2D:
%   Retorna o valor do pseudo-espectro MUSIC (8x8) em (xx,yy)
%   -> Aqui usamos a forma 1 / |a'*(Un*Un')*a|, mas
%      no caso abaixo, vamos "maximizar" => definimos nm_val = -valPS
%      (ou seja, "ps_func" = "valPS").

[Nx,Ny] = size(elements_x);
N = Nx*Ny;
a_ = zeros(N,1);

d_ref_ = sqrt( (x_ref - xx)^2 + (y_ref - yy)^2 + z_ref^2 );
idx=1;
for ix_=1:Nx
    for iy_=1:Ny
        dx_ = elements_x(ix_, iy_);
        dy_ = elements_y(ix_, iy_);
        dz_ = elements_z(ix_, iy_);
        d_km = sqrt( (dx_-xx)^2 + (dy_-yy)^2 + dz_^2 );
        ph_ = -(2*pi/lambda)*( d_ref_ - d_km );
        a_(idx) = exp(1j*ph_);
        idx=idx+1;
    end
end

% Valor no denominador do pseudo-espectro:  abs(a^H * (Un*Un') * a).
den_ = abs(a_'*(Un_2D*Un_2D')*a_);
if den_<1e-12
    valPS = 1e6;  % para evitar singularidade
else
    valPS = 1 / den_;
end
end

function [Jxx, Jyy, Jxy] = precomputeCRBterms(elements_x, elements_y, ...
    elements_z, user_xyz, lambda)
% Cálculo do CRB (soma das derivadas parciais) usando eq. teórica
user_x = user_xyz(1);
user_y = user_xyz(2);
user_z = user_xyz(3);

Nant = numel(elements_x);
Jxx=0; Jyy=0; Jxy=0;

x_ref = elements_x(1);
y_ref = elements_y(1);
z_ref = elements_z(1);
d_ref = sqrt( (x_ref-user_x)^2 + (y_ref-user_y)^2 + (z_ref-user_z)^2 );

for i=1:Nant
    x_i  = elements_x(i);
    y_i  = elements_y(i);
    z_i  = elements_z(i);
    d_i  = sqrt( (x_i-user_x)^2 + (y_i-user_y)^2 + (z_i-user_z)^2 );

    delta_i   = (user_x - x_i);
    delta_ref = (user_x - x_ref);

    partial_x_1 = (delta_i / d_i^3);
    partial_x_2 = (2*pi / lambda)^2 * ( (delta_i)/(d_i^2) - (delta_ref)/(d_ref*d_i) )^2;
    px = partial_x_1^2 + partial_x_2;  
    Jxx = Jxx + px;

    partial_y_1 = (user_y)/(d_i^3);
    partial_y_2 = (2*pi/lambda)^2 * ( (user_y)/(d_i^2) - (user_y)/(d_ref*d_i) )^2;
    py = partial_y_1^2 + partial_y_2;
    Jyy = Jyy + py;

    cross_part_1 = (delta_i / d_i^3)*(user_y / d_i^3);
    cross_part_2 = (2*pi/lambda)^2 * ...
        ( (delta_i)/(d_i^2) - (delta_ref)/(d_ref*d_i) ) * ...
        ( (user_y)/(d_i^2) - (user_y)/(d_ref*d_i) );
    Jxy = Jxy + (cross_part_1 + cross_part_2);
end
end
