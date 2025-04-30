function cramer_rao_bound = crb(L, URA, UEs, lambda, noise_power)
    % Fator de normalizacao da FIM (conforme CRB derivado)
    P_norm = L * (lambda / (sqrt(2)*pi*sqrt(noise_power)))^2;
    
    % Calculo dos termos da FIM nao escalados 
    [Jxx_c, Jyy_c, Jxy_c] = precomputeCRBterms(URA, UEs, lambda);
    
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
    cramer_rao_bound = sqrt(crb_x + crb_y);
end

function [Jxx, Jyy, Jxy] = precomputeCRBterms(URA, user_xyz, lambda)
% Calcula os termos não escalados da matriz de informação de Fisher (FIM)
% URA: matriz M×3 com coordenadas [x, y, z]
% user_xyz: vetor [x_u, y_u, z_u]

    user_x = user_xyz(1);
    user_y = user_xyz(2);
    user_z = user_xyz(3);

    x_ref = URA(1,1);
    y_ref = URA(1,2);
    z_ref = URA(1,3);

    d_ref = sqrt((x_ref - user_x)^2 + (y_ref - user_y)^2 + (z_ref - user_z)^2);

    Jxx = 0; Jyy = 0; Jxy = 0;

    for i = 1:size(URA,1)
        x_i = URA(i,1);
        y_i = URA(i,2);
        z_i = URA(i,3);

        d_i = sqrt((x_i - user_x)^2 + (y_i - user_y)^2 + (z_i - user_z)^2);

        delta_i   = user_x - x_i;
        delta_ref = user_x - x_ref;

        % Derivadas em x
        partial_x_1 = delta_i / d_i^3;
        partial_x_2 = (2*pi / lambda)^2 * ((delta_i/d_i^2) - (delta_ref/(d_ref*d_i)))^2;
        px = partial_x_1^2 + partial_x_2;
        Jxx = Jxx + px;

        % Derivadas em y
        partial_y_1 = user_y / d_i^3;
        partial_y_2 = (2*pi / lambda)^2 * ((user_y/d_i^2) - (user_y/(d_ref*d_i)))^2;
        py = partial_y_1^2 + partial_y_2;
        Jyy = Jyy + py;

        % Termo cruzado x-y
        cross_part_1 = (delta_i / d_i^3) * (user_y / d_i^3);
        cross_part_2 = (2*pi / lambda)^2 * ...
            ((delta_i/d_i^2) - (delta_ref/(d_ref*d_i))) * ...
            ((user_y/d_i^2) - (user_y/(d_ref*d_i)));
        Jxy = Jxy + (cross_part_1 + cross_part_2);
    end
end