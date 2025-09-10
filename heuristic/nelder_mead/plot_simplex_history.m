function plot_simplex_history(Pmusic_fixo, x_grid, y_grid, simplex_history)
    [X, Y] = meshgrid(x_grid, y_grid);
    Z = Pmusic_fixo';

    figure;
    contourf(X, Y, Z, 40); colorbar;
    hold on;

    for k = 1:length(simplex_history)
        simplex = simplex_history{k};
        fill(simplex(:,1), simplex(:,2), 'r', 'FaceAlpha', 0.2, ...
            'EdgeColor', 'k', 'LineWidth', 1.2);
        plot(mean(simplex(:,1)), mean(simplex(:,2)), 'ko', 'MarkerSize', 6, ...
            'MarkerFaceColor', 'w');
        pause(0.3);
    end

    title('Evolucao do Simplex no Pseudoespectro MUSIC');
    xlabel('x (m)');
    ylabel('y (m)');
end