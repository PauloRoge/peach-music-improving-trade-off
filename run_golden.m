startup;
tic;
% Geracao da geometria das ULAS
[ULA_h, ULA_v, x_h, x_v, z_h, z_v] = subarrays_ula(8, 8, d_x, d_z, elev, lambda);
ref = ULA_h(1,:);  % ponto de referência

% % SINAL RECEBIDO (modelo fisico fiel) baseado na SNR
[Yh, Yv, Y] = signals_ula(pos, ULA_h, ULA_v, lambda, L, alpha, SNR_dB, P_tx);

%%%%%%%%%%%%%%%%%%%%%%%%% PEAK FINDER %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Funcao anonima de resposta de array
responsearray = @(x, y, z) steering_vector(Mx, Mz, elev, d_x, d_z, lambda, x, y, z);
[Pmusic] = pseudospectrum(responsearray, Y, L);

[Un_h, Un_v, pos_est] = peach_aurea( ...
        Yh, Yv, L, x, n_hiper, x_h, z_h, x_v, z_v, ref, lambda, y, n_circ, pos);

% Exibicao dos resultados
fprintf('\nPosicao PEACH-MUSIC: (%.2f, %.2f)', pos_est(1), pos_est(2));

% Calculo dos erros euclidianos entre estimativa e posicao real
erro_peach  = norm(pos_est - pos(1:2));

% Exibicao formatada dos resultados
fprintf('\n=============================================================\n');
fprintf('Erro Euclidiano Golden SectionPEACH        = %.3f m\n', erro_peach);
fprintf('===============================================================\n');

% Exibir resultado
fprintf('Posicao REAL do usuario: [%.2f, %.2f]\n', pos(1), pos(2));
fprintf('Posicao PEACH estimada: [%.2f, %.2f]\n', pos_est(1), pos_est(2));
toc;

% crb_eucl = crb_vetorizado(L, URA, UEs, lambda, P_tx, SNR_dB);
% fprintf('\nCRB (Euclidiano) = %.4f m\n', crb_eucl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       FUNCOES LOCAIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%-------------------------------------------------------------------------
% FUNÇÃO PARA CALCULAR O PSEUDO-ESPECTRO
%-------------------------------------------------------------------------
function value = pmusic(responsearray, pos, Vn)
    % Compatibiliza entrada 2D ou 3D

    a = responsearray(pos(1), pos(2), pos(3));
    value = 1 / abs(a' * (Vn * Vn') * a);
end