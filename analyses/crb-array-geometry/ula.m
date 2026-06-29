function [ULA] = ula(M, d, elev, lambda, plotflag)
%   [ULA] = subarray_ula(M, d, elev, lambda, plotflag)
%
% Parâmetros:
%   M         : número de elementos da ULA
%   d         : espaçamento entre elementos (m)
%   elev      : elevação (offset em z)
%   lambda    : comprimento de onda (para escala no gráfico)
%   plotflag  : true/false para plotar (default = false)
%
% Saída:
%   ULA : [M x 3] posições (x, y, z), com y = 0

    if nargin < 5
        plotflag = false;
    end

    % Gera coordenadas igualmente espaçadas em torno do centro
    x = ((0:M-1) - (M-1)/2) * d;
    z = elev * ones(1, M);
    y = zeros(1, M);

    % Matriz de posições
    ULA = [x(:), y(:), z(:)];

    % Plotagem opcional 
    if plotflag
        figure; hold on; grid on;
        plot3(ULA(:,1), ULA(:,2), ULA(:,3), 'o-', ...
              'MarkerFaceColor', [0.1, 0.4, 0.8], 'MarkerEdgeColor', 'k', ...
              'LineWidth', 1.5, 'DisplayName', 'ULA');
        xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
        title('Posição dos Elementos da ULA');
        legend('Location','best'); axis equal; view(45,25);
        xlim([-(M/2)*d - lambda, (M/2)*d + lambda]);
        ylim([-0.05, 0.05]);
        zlim([elev - lambda, elev + lambda]);
        set(gca, 'FontSize', 12);
    end
end
