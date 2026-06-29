function cramer_rao_bound = crb_vetorizado(L, URA, pos, lambda, P_tx, SNR_dB)
    % Conversão de SNR para escala linear
    SNR_lin = 10^(SNR_dB/10);

    % -------------------------
    % Potência total recebida
    % -------------------------
    diff = URA - pos;
    d_m = sqrt(sum(diff.^2, 2));  % vetor de distâncias
    beta = (lambda ./ (4*pi*d_m)).^2;
    total_rx_power = P_tx * sum(beta);

    % Potência do ruído
    noise_power = total_rx_power / SNR_lin;

    % Fator de normalização do FIM
    P_norm = L * (lambda / (sqrt(2)*pi*sqrt(noise_power)))^2;

    % -------------------------
    % Termos da FIM
    % -------------------------
    [Jxx_c, Jyy_c, Jxy_c] = CRBterms_vetorizado(URA, pos, lambda);
    
    Jxx_s = P_norm * Jxx_c;
    Jyy_s = P_norm * Jyy_c;
    Jxy_s = P_norm * Jxy_c;

    % Denominadores da inversa da FIM
    den_x = Jxx_s - (Jxy_s^2 / Jyy_s);
    den_y = Jyy_s - (Jxy_s^2 / Jxx_s);

    % CRBs individuais
    crb_x = 1 / den_x;
    crb_y = 1 / den_y;

    % CRB combinado (erro euclidiano mínimo)
    cramer_rao_bound = sqrt(crb_x + crb_y);
end

function [Jxx, Jyy, Jxy] = CRBterms_vetorizado(URA, user_xyz, lambda)
    % Extração dos dados
    x_i = URA(:,1);
    y_i = URA(:,2);
    z_i = URA(:,3);
    
    x_u = user_xyz(1);
    y_u = user_xyz(2);
    z_u = user_xyz(3);

    % Referência (primeira antena)
    x_ref = x_i(1);
    y_ref = y_i(1);
    z_ref = z_i(1);

    % Distâncias
    d_i = sqrt((x_i - x_u).^2 + (y_i - y_u).^2 + (z_i - z_u).^2);
    d_ref = sqrt((x_ref - x_u)^2 + (y_ref - y_u)^2 + (z_ref - z_u)^2);

    % Deltas
    delta_i = x_u - x_i;
    delta_ref = x_u - x_ref;

    % Termos Jxx
    partial_x_1 = delta_i ./ d_i.^3;
    partial_x_2 = (2*pi/lambda)^2 .* ((delta_i ./ d_i.^2 - delta_ref ./ (d_ref .* d_i)).^2);
    px = partial_x_1.^2 + partial_x_2;
    Jxx = sum(px);

    % Termos Jyy
    partial_y_1 = y_u ./ d_i.^3;
    partial_y_2 = (2*pi/lambda)^2 .* ((y_u ./ d_i.^2 - y_u ./ (d_ref .* d_i)).^2);
    py = partial_y_1.^2 + partial_y_2;
    Jyy = sum(py);

    % Termos Jxy
    cross_part_1 = (delta_i .* y_u) ./ d_i.^6;
    cross_part_2 = (2*pi/lambda)^2 .* ...
                   ((delta_i ./ d_i.^2 - delta_ref ./ (d_ref .* d_i)) .* ...
                    (y_u ./ d_i.^2 - y_u ./ (d_ref .* d_i)));
    Jxy = sum(cross_part_1 + cross_part_2);
end
