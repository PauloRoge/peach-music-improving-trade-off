close all;
clc;
startup;
tic;
% Geracao da geometria da URA
[URA, URA_x, URA_z, x_h, x_v, z_h, z_v] = subarrays(Mx, Mz, ...
    d_x, d_z, elev, lambda, plt_array);

URA_y = zeros(size(URA_x));

% % SINAL RECEBIDO (modelo fisico fiel) baseado na SNR
[Yh, Yv, Y] = signals_los(UEs, URA, lambda, L, alpha, ...
    SNR_dB, P_tx, Mx, Mz);

%%%%%%%%%%%%%%%%%%%%%%%%% PEAK FINDER %%%%%%%%%%%%%%%%%%%%%%%%%%
% Elemento de referência
ref = URA(1,:);

% Funcao anonima de resposta de array
responsearray = @(x, y, z) steering_vector(Mx, Mz, elev, ...
    d_x, d_z, lambda, x, y, z);

[Pmusic] = pseudospectrum(responsearray, Y, L);

%Estimativa PEACH com subarranjos
[Un_h, Un_v, est_peach] = peach_golden(Yh, Yv, L, x, n_hiper, ...
    x_h, z_h, x_v, z_v, ref, lambda, y, n_circ, pos);

% Exibicao dos resultados
fprintf('\nPosicao PEACH-MUSIC: (%.2f, %.2f)', ...
    est_peach(1), est_peach(2));
 
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

[nm_est, history] = nelder_mead(URA, est_peach, Un, ...
    lambda, ref, deltaArea, numIterNM, 1e-6, true, x, y);

% Calculo dos erros euclidianos entre estimativa e posicao real
erro_peach  = norm(est_peach - pos(1:2));
erro_nelder = norm(nm_est - pos(1:2));

printResults;

crb_eucl = crb_vetorizado(L, URA, UEs, lambda, P_tx, SNR_dB);

fprintf('\nCRB (Euclidiano) = %.4f m\n', crb_eucl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       FUNCOES LOCAIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------
% FUNCTION STEERING VECTOR
% -------------------------------------------------------------
function a = steering_vector(Mx, Mz, elev, d_x, d_z, lambda, ...
    x, y, z)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Pmusic, Un] = pseudospectrum(responsearray, Y, ...
    snapshots)
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

%------------------------------------------------------------
% FUNÇÃO PARA CALCULAR O PSEUDO-ESPECTRO
%------------------------------------------------------------
function value = pmusic(responsearray, pos, Vn)
    % Compatibiliza entrada 2D ou 3D

    a = responsearray(pos(1), pos(2), pos(3));
    value = 1 / abs(a' * (Vn * Vn') * a);
end