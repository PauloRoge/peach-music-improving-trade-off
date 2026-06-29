function est = music(Y, URA, L, lambda, ref, x_grid, y_grid)
% MUSIC   Estima 2D [(x,y)] via MUSIC usando URA completa
%
%   est = music(Y, URA, L, lambda, ref, x_grid, y_grid) retorna um vetor [x_est, y_est]
%   que maximiza o pseudospectro de MUSIC. 
%   - Y       : M x L matriz de snapshots (M = Mx*Mz elementos da URA)
%   - URA     : M x 3, cada linha URA(k,:) = [x_k, y_k, z_k]
%   - L       : número de snapshots
%   - lambda  : comprimento de onda
%   - ref     : 1 x 3, [x_ref, y_ref, z_ref] coordenadas do elemento referência
%   - x_grid  : vetor de candidatos em x (1 x n_x)
%   - y_grid  : vetor de candidatos em y (1 x n_y)
%
%   Retorna:
%   - est     : [x_est, y_est], estimativa 2D de posição do usuário 
%

    % Número de elementos M
    [M, ~] = size(URA);
    
    % 1) Matriz de covariância
    R = (Y * Y') / L;
    
    % 2) Autovalores e autovetores
    [V, D] = eig(R);
    [~, ind] = sort(diag(D), 'descend');
    V = V(:, ind);
    
    % 3) Assumindo 1 sinal, há (M-1) autovetores no subespaço-ruído
    Un = V(:, 2:end);   % tamanho M x (M-1)
    
    % 4) Pré-alocar pisc
    n_x = length(x_grid);
    n_y = length(y_grid);
    Pmusic = zeros(n_x, n_y);
    
    % Coordenadas de referência
    x_ref = ref(1);
    y_ref = ref(2);
    z_ref = ref(3);
    
    % 5) Para cada ponto candidato (x_i, y_j), calcule steering vector
    for ix = 1:n_x
        x_i = x_grid(ix);
        for jy = 1:n_y
            y_j = y_grid(jy);
            
            % Supondo usuário em z0 = 0; se for outra cota, basta ajustar aqui
            z0 = 0;
            
            % Distância do ponto candidato ao elemento de referência
            d_ref = sqrt( (x_i - x_ref)^2 + (y_j - y_ref)^2 + (z0 - z_ref)^2 );
            
            % Montar vetor de resposta a
            a = zeros(M,1);
            for k = 1:M
                xk = URA(k,1);
                yk = URA(k,2);
                zk = URA(k,3);
                d_k = sqrt( (x_i - xk)^2 + (y_j - yk)^2 + (z0 - zk)^2 );
                phi = - (2*pi / lambda) * (d_ref - d_k);
                a(k) = exp(1j * phi);
            end
            
            % 6) Avaliar pseudospectro de MUSIC em (x_i, y_j)
            denom = a' * (Un * Un') * a;
            Pmusic(ix,jy) = 1 / abs(denom);
        end
    end
    
    % 7) Encontrar índice do pico máximo
    [~, idx_max] = max(Pmusic(:));
    [ix_max, jy_max] = ind2sub(size(Pmusic), idx_max);
    
    % 8) Coordenadas estimadas
    est = [x_grid(ix_max), y_grid(jy_max)];
end
