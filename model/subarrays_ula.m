function [ULA_h, ULA_v, x_h, x_v, z_h, z_v] = subarrays_ula(Mx, Mz, d_x, d_z, elev, lambda)
    % ------------------- geometria -------------------
    x_offset = (Mx-1)*d_x/2;          % centralizar array em x
    x_all    = (0:Mx-1).' * d_x - x_offset;   % Mx×1
    z_all    = (0:Mz-1)   * d_z + elev;       % 1×Mz

    % ULA horizontal (linha 1 da malha)
    x_h = x_all;                      
    z_h = elev * ones(Mx,1);          
    ULA_h = [x_h, zeros(Mx,1), z_h];

    % ULA vertical (coluna 1 da malha)
    x_v = -x_offset * ones(Mz,1);     
    z_v = z_all.';                     
    ULA_v = [x_v, zeros(Mz,1), z_v];

    %chama o plot passando os vetores já prontos
    %plot_subarrays_ula(ULA_h, ULA_v, Mx, Mz, d_x, d_z, elev, lambda);
end

% function plot_subarrays_ula(ULA_h, ULA_v, Mx, Mz, d_x, d_z, elev, lambda)
%     figure; hold on; grid on;
% 
%     % ULA horizontal
%     plot3(ULA_h(:,1), ULA_h(:,2), ULA_h(:,3), 's', ...
%           'MarkerEdgeColor','k','MarkerFaceColor',[0,0.2,0.6], ...
%           'DisplayName','ULA Horizontal');
% 
%     % ULA vertical
%     plot3(ULA_v(:,1), ULA_v(:,2), ULA_v(:,3), 'o', ...
%           'MarkerEdgeColor','k','MarkerFaceColor',[0.6,0.2,0], ...
%           'DisplayName','ULA Vertical');
% 
%     % Configurações do gráfico
%     xlabel('x [m]');
%     ylabel('y [m]');
%     zlabel('z [m]');
%     legend('Location','best');
%     title('Posição das Antenas no Espaço (Subarrays)');
%     view(45, 25);
%     axis equal;
% 
%     % Ajusta os limites com margem
%     xlim([-(Mx/2)*d_x - lambda, (Mx/2)*d_x + lambda]);
%     ylim([-0.05, 0.05]);             
%     zlim([elev - lambda, elev + (Mz+1)*d_z]);
%     set(gca, 'FontSize', 12);
% end
