function [URA, URA_x, URA_z, x_h, x_v, z_h, z_v] = subarrays(Mx, Mz, d_x, d_z, ...
elev, lambda, plot)
    x_offset = (Mx-1)*d_x/2;

    elements_x_a = zeros(Mx, Mz);
    elements_z_a = zeros(Mx, Mz);
    
    for i = 1:Mx
        for j = 1:Mz
            elements_x_a(i,j) = (i-1)*d_x - x_offset;
            elements_z_a(i,j) = (j-1)*d_z + elev;
        end
    end

    % Linha fixa (j = 1) -> ULA horizontal
    x_h = elements_x_a(:,1);   % Varia em i
    z_h = elements_z_a(:,1);

    % Coluna fixa (i = 1) -> ULA vertical
    x_v = elements_x_a(1,:)';  % Varia em j
    z_v = elements_z_a(1,:)';
    
    URA_x = elements_x_a;
    URA_z = elements_z_a;

    % vetorizar mudando de linha apos varrer todas as colunas
    URA_x_vetorizado = reshape(elements_x_a.', [], 1); 
    URA_z_vetorizado = reshape(elements_z_a.', [], 1);
    URA_y_vetorizado = zeros(size(URA_x_vetorizado));
    URA = [URA_x_vetorizado, URA_y_vetorizado, URA_z_vetorizado];

    if (plot == true)
        disp(['URA size: ', num2str(size(URA_x, 1)), ' x ', ...
            num2str(size(URA_x, 2))])
        disp(['Mx length: ', num2str(length(x_h))])
        disp(['Mz length: ', num2str(length(x_v))])

        % Prepara o grafico
        figure;
        hold on;
        grid on;

        % Plota todas as antenas do URA (preto, contorno)
        plot3(URA_x(:), zeros(size(URA_x(:))), URA_z(:), 'ko', ...
            'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'URA');

        % ULA vertical (Yh) – marcador circular azul preenchido
        plot3(x_h, zeros(size(x_h)), z_h, 's', ...
            'MarkerEdgeColor', 'black', 'MarkerFaceColor', ...
            [0.0, 0.2, 0.6], 'MarkerSize', 12, 'LineWidth', 1.5, ...
            'DisplayName', 'ULA 1');

        % ULA horizontal (Yv) – marcador circular amarelo preenchido
        plot3(x_v, zeros(size(x_v)), z_v, 's', ...
            'MarkerEdgeColor', 'black', 'MarkerFaceColor', ...
            [0.0, 0.2, 0.6], 'MarkerSize', 12, 'LineWidth', 1.5, ...
            'DisplayName', 'ULA 2');

        % ---------------------------------------------------------------
        % Intersecao entre subarrays (ponto comum Yh ? Yv) ? marcador verde
        % ---------------------------------------------------------------
        % Normaliza vetores em coluna
        xh_col = [x_h(:), z_h(:)];
        xv_col = [x_v(:), z_v(:)];

        % Intersecao com tolerancia numerica
        intersec = [];
        for i = 1:size(xh_col,1)
            for j = 1:size(xv_col,1)
                if norm(xh_col(i,:) - xv_col(j,:)) < 1e-6
                    intersec = xh_col(i,:);
                    break;
                end
            end
        end

        % Plota intersecao se existir
        if ~isempty(intersec)
            plot3(intersec(1), 0, intersec(2), 'o', ...
                'MarkerEdgeColor', 'black', 'MarkerFaceColor', [1.0, 0.5, 0.0], ...
                'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'Ref');
        end

        % Configuracoes do grafico
        xlabel('x [m]');
        ylabel('y [m]');
        zlabel('z [m]');
        legend('Location', 'best');
        title('Posicao das Antenas no Espaco (URA + Subarrays)');
        view(45, 25);
        axis equal;
        
        % Ajusta os limites com margem
        xlim([-(Mx/2)*d_x - lambda, (Mx/2)*d_x + lambda]);
        ylim([-0.05, 0.05]);  % y fixo para visualizacao 3D
        zlim([elev - lambda, elev + (Mz+1)*d_z]);
        set(gca, 'FontSize', 12);

    end
end