tic;

% Geracao da geometria da URA
[URA, URA_x, URA_z, x_h, x_v, z_h, z_v] = subarrays(Mx, Mz, d_x, d_z, elev, lambda, plt_array);
URA_y = zeros(size(URA_x));  % plano XZ ? y fixo = 0 para todas as antenas

% % SINAL RECEBIDO (modelo fisico fiel) baseado na pot.tx e no ruído direto 
%  [Yh, Yv, Y] = signals(UEs, URA_x, URA_y, URA_z, lambda, L, alpha, ...
%      P_tx, N_dBm);

% % SINAL RECEBIDO (modelo fisico fiel) baseado na SNR
[Yh, Yv, Y] = signals(UEs, URA, lambda, L, alpha, SNR_dB, P_tx, Mx, Mz);
%%%%%%%%%%%%%%%%%%%%%%%%% PEAK FINDER %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Elemento de referência
ref = URA(1,:);

% Funcao anonima de resposta de array
responsearray = @(x, y, z) steering_vector(Mx, Mz, elev, d_x, d_z, lambda, x, y, z);
[Pmusic] = pseudospectrum(responsearray, Y, L);

% Estimativa PEACH com subarranjos
[Un_h, Un_v, pos_est] = peach_analitico(Yh, Yv, L, x, n_hiper, ...
    x_h, z_h, x_v, z_v, ref, ...
    lambda, y, n_circ, pos);

% Exibicao dos resultados
fprintf('\nPosicao PEACH-MUSIC: (%.2f, %.2f)', pos_est(1), pos_est(2));
 
% MUSIC (completo)
 %-----------------------------------------------
    % DIVISÃO DO SUBESPAÇO (Vn)
    %-----------------------------------------------
    Cov = (Y * Y') / L;
    [eigenvectors, eigenvalues] = eig(Cov); 
    estimated_sources = 1;
    [~, i] = sort(diag(eigenvalues), 'descend'); 
    eigenvectors = eigenvectors(:, i);
    Un = eigenvectors(:, estimated_sources+1:end);
    %-----------------------------------------------

% Parâmetros para refinamento Nelder-Mead
deltaArea = 5;
numIterNM = 100;

% [nm_est, simplex_history] = nelder_mead(URA, pos_est, Un, lambda, ref, ...
%     deltaArea, numIterNM, x, y);

[nm_est, simplex_history] = nelder_mead(URA, pos_est, Un, lambda, ref, ...
    deltaArea, numIterNM, 1e-6, true, x, y);


% Calculo dos erros euclidianos entre estimativa e posicao real
erro_peach  = norm(pos_est - pos(1:2));
erro_nelder = norm(nm_est - pos(1:2));

% Exibicao formatada dos resultados
fprintf('\n====================================================\n');
fprintf('Erro Euclidiano PEACH        = %.3f m\n', erro_peach);
fprintf('Erro Euclidiano Nelder-Mead  = %.7f m\n', erro_nelder);
fprintf('====================================================\n');

% Exibir resultado
fprintf('Posicao REAL do usuario:     [%.2f, %.2f]\n', pos(1), pos(2));
fprintf('Posicao PEACH estimada:      [%.2f, %.2f]\n', pos_est(1), pos_est(2));
fprintf('Posicao Nelder-Mead refinada:[%.2f, %.2f]\n', nm_est(1), nm_est(2));
toc;

if (plt_spectrum == 1)
    % Movimentacao do Simplex
    % Gera pseudoespectro fixo (para visualizacao)
    Pmusic_fixo = zeros(length(x_grid), length(y_grid));
    for i = 1:length(x_grid)
        for j = 1:length(y_grid)
            pos = [x_grid(i), y_grid(j), 0];
            a = responsearray(pos(1), pos(2), 0);
            Pmusic_fixo(i,j) = Pmusic(pos);
        end
    end

    % Normaliza
    Pmusic_fixo = abs(Pmusic_fixo);
    Pmusic_fixo = Pmusic_fixo ./ max(Pmusic_fixo(:));

    figure;
    imagesc(x_grid, y_grid, 10*log10(Pmusic_fixo'));
    axis xy;
    xlabel('x (m)');
    ylabel('y (m)');
    title(sprintf('Pseudoespectro MUSIC (SNR = %d dB), %d×%d URA', SNR_dB, Mx, Mz));
    colorbar;

    % [X, Y] = meshgrid(x_grid, y_grid);
    % figure;
    % surf(X, Y, 10*log10(Pmusic_fixo'));
    % xlabel('x (m)');
    % ylabel('y (m)');
    % zlabel('Pseudoespectro (dB)');
    % title('Pseudoespectro MUSIC');
    % shading interp;
    % colorbar;
    % view(45, 45); % Ângulo de visão


    if(plt_neldermead == 1)
    % Visualizacao
    plot_simplex_history(Pmusic_fixo, x_grid, y_grid, simplex_history);
    end
end

x_ref = URA_x(1,1);
y_ref = URA_y(1,1);
z_ref = URA_z(1,1);

x_lim = 45;
y_lim = 45;
y_min = 10;

% Defina a resolução e a faixa de varredura:
Nx_plot = 100;  
Ny_plot = 100;

% Vetores de busca no plano (x, y)
x_vec = linspace(-x_lim, x_lim, Nx_plot);
y_vec = linspace(0, y_lim, Ny_plot);

% Inicialização dos mapas
Pmap_h = zeros(Ny_plot, Nx_plot);
Pmap_v = zeros(Ny_plot, Nx_plot);

% Índices dos subarranjos
indices_h = ((1:Mx) - 1)*Mz + 1;   % Linha fixa (coluna j=1)
indices_v = 1:Mz;                 % Coluna fixa (linha i=1)

% Extração das posições dos subarrays
URA_h = URA(indices_h, :);   % Subarray horizontal (Mx × 3)
URA_v = URA(indices_v, :);   % Subarray vertical   (Mz × 3)

% --- Mapa de Calor: Subarranjo Horizontal ---
for i = 1:Nx_plot
    for j = 1:Ny_plot
        x_cand = x_vec(i);
        y_cand = y_vec(j);

        % Distância de referência (antena 1)
        d_ref = sqrt((ref(1) - x_cand)^2 + (ref(2) - y_cand)^2 + ref(3)^2);

        % Distâncias dos elementos do subarray
        d_km = sqrt((URA_h(:,1) - x_cand).^2 + ...
                    (URA_h(:,2) - y_cand).^2 + URA_h(:,3).^2);

        % Vetor de steering
        phase = -(2*pi/lambda) * (d_ref - d_km);
        ah = exp(1j * phase);

        % Pseudoespectro
        Pmap_h(j, i) = 1 / abs(ah' * (Un_h * Un_h') * ah);
    end
end

% --- Mapa de Calor: Subarranjo Vertical ---
for i = 1:Nx_plot
    for j = 1:Ny_plot
        x_cand = x_vec(i);
        y_cand = y_vec(j);

        d_ref = sqrt((ref(1) - x_cand)^2 + (ref(2) - y_cand)^2 + ref(3)^2);

        d_km = sqrt((URA_v(:,1) - x_cand).^2 + ...
                    (URA_v(:,2) - y_cand).^2 + URA_v(:,3).^2);

        phase = -(2*pi/lambda) * (d_ref - d_km);
        a_v = exp(1j * phase);

        Pmap_v(j, i) = 1 / abs(a_v' * (Un_v * Un_v') * a_v);
    end
end


Pmap_sum = Pmap_h + Pmap_v;
Pmap_mul = Pmap_h .* Pmap_v;

if(plt_peach)
    figure('Name','Mapas de Calor - Somente y >= 0','Color',[1 1 1]);
    
    % Subplot 1
    subplot(2,2,1);
    imagesc(x_vec, y_vec, 10*log10(Pmap_h)/ max(max( 10*log10(Pmap_h))));
    set(gca,'YDir','normal'); colorbar;
    title('Horizontal (Hipérbole)');
    xlabel('x (m)'); ylabel('y (m)');
    hold on; plot(UEs(1), UEs(2), 'ro','MarkerFaceColor','r','MarkerSize',7);
    hold off;
    
    % Subplot 2
    subplot(2,2,2);
    imagesc(x_vec, y_vec, 10*log10(Pmap_v)/max(max( 10*log10(Pmap_v))));
    set(gca,'YDir','normal'); colorbar;
    title('Vertical (Círculo)');
    xlabel('x (m)'); ylabel('y (m)');
    hold on; plot(UEs(1), UEs(2), 'ro','MarkerFaceColor','r','MarkerSize',7);
    hold off;
    
    % Subplot 3
    subplot(2,2,3);
    imagesc(x_vec, y_vec, 10*log10(Pmap_sum)/max(max( 10*log10(Pmap_sum))));
    set(gca,'YDir','normal'); colorbar;
    title('Soma (Hipérbole + Círculo)');
    xlabel('x (m)'); ylabel('y (m)');
    hold on; plot(UEs(1), UEs(2), 'ro','MarkerFaceColor','r','MarkerSize',7);
    hold off;
    
    % Subplot 3
    subplot(2,2,4);
    imagesc(x_vec, y_vec, 10*log10(Pmap_mul)/max(max( 10*log10(Pmap_mul))));
    set(gca,'YDir','normal'); colorbar;
    title('Soma (Hipérbole + Círculo)');
    xlabel('x (m)'); ylabel('y (m)');
    hold on; plot(UEs(1), UEs(2), 'ro','MarkerFaceColor','r','MarkerSize',7);
    hold off;
end

if(plt_hiper)
    % Subplot 1
    figure(1337)
    imagesc(x_vec, y_vec, 10*log10(Pmap_h)/ max(max( 10*log10(Pmap_h))));
    set(gca,'YDir','normal'); colorbar;
    xlabel('x (m)'); ylabel('y (m)');
    hold on; plot(UEs(1), UEs(2), 'ro','MarkerFaceColor','r','MarkerSize',7);
    hold off;
end

if(plt_circle)
% Subplot 2
    figure(1338)
    imagesc(x_vec, y_vec, 10*log10(Pmap_v)/max(max( 10*log10(Pmap_v))));
    set(gca,'YDir','normal'); colorbar;
    xlabel('x (m)'); ylabel('y (m)');
    hold on; plot(UEs(1), UEs(2), 'ro','MarkerFaceColor','r','MarkerSize',7);
    hold off;
end

if(plt_itersec)
    % Subplot 3
    figure(1339)
    imagesc(x_vec, y_vec, 10*log10(Pmap_sum)/max(max( 10*log10(Pmap_sum))));
    set(gca,'YDir','normal'); colorbar;
    xlabel('x (m)'); ylabel('y (m)');
    hold on; plot(UEs(1), UEs(2), 'ro','MarkerFaceColor','r','MarkerSize',7);
    hold off;
end

% Calculo do CRB teorico (normalizado) com URA direto
% Calculo dos termos da FIM nao escalados 
[Jxx_c, Jyy_c, Jxy_c] = precomputeCRBterms(URA, UEs, lambda);

% Calculo da potencia recebida total (usado para noise power)
total_rx_power = 0;
for i = 1:size(URA, 1)
    d_i = norm(URA(i,:) - UEs);  % distância do usuário até a antena i
    beta_i = (lambda / (4*pi*d_i))^2;
    total_rx_power = total_rx_power + power * beta_i;
end

% Potencia do ruído
SNR_lin = 10^(SNR_dB/10);
noise_power = total_rx_power / SNR_lin;

% Fator de normalizacao da FIM (conforme CRB derivado)
P_norm = L * (lambda / (sqrt(2)*pi*sqrt(noise_power)))^2;

% Termos escalados
Jxx_s = P_norm * Jxx_c;
Jyy_s = P_norm * Jyy_c;
Jxy_s = P_norm * Jxy_c;

% Denominadores do CRB
den_x = Jxx_s - (Jxy_s^2 / Jyy_s);
den_y = Jyy_s - (Jxy_s^2 / Jxx_s);

% CRB individuais
crb_x = 1 / den_x;
crb_y = 1 / den_y;

% CRB combinado euclidiano
crb_eucl = sqrt(crb_x + crb_y);

fprintf('\nCRB (Euclidiano) = %.4f m\n', crb_eucl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       FUNCOES LOCAIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------------------------------------------------
% FUNCTION STEERING VECTOR
% ------------------------------------------------------------------
function a = steering_vector(Mx, Mz, elev, d_x, d_z, lambda, x, y, z)
    x_pos = (0:Mx-1) * d_x;
    z_pos = (0:Mz-1) * d_z + elev;
    
    a = zeros(Mx * Mz, 1);
    index = 1;

    % Posicao do elemento de referencia (1,1)
    x_ref = x_pos(1);
    z_ref = z_pos(1);
    d_ref = sqrt((x - x_ref)^2 + y^2 + (z - z_ref)^2);

    for i = 1:Mx
        for j = 1:Mz
            x_k = x_pos(i);
            z_k = z_pos(j);
            d_k = sqrt((x - x_k)^2 + y^2 + (z - z_k)^2);
            phase_k = -(2*pi/lambda)*(d_ref - d_k);
            a(index) = exp(1j * phase_k);
            index = index + 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Pmusic, Un] = pseudospectrum(responsearray, Y, snapshots)
    %-----------------------------------------------
    % DIVISÃO DO SUBESPAÇO (Vn)
    %-----------------------------------------------
    Cov = (Y * Y') / snapshots;
    [eigenvectors, eigenvalues] = eig(Cov); 
    estimated_sources = 1;
    [~, i] = sort(diag(eigenvalues), 'descend'); 
    eigenvectors = eigenvectors(:, i);
    Vn = eigenvectors(:, estimated_sources+1:end);
    %-----------------------------------------------
    Un = Vn;
    Pmusic = @(pos) pmusic(responsearray, pos, Vn);
end

%-------------------------------------------------------------------------------------------------------
% FUNÇÃO PARA CALCULAR O PSEUDO-ESPECTRO
%-------------------------------------------------------------------------------------------------------
function value = pmusic(responsearray, pos, Vn)
    % Compatibiliza entrada 2D ou 3D

    a = responsearray(pos(1), pos(2), pos(3));
    value = 1 / abs(a' * (Vn * Vn') * a);
end



