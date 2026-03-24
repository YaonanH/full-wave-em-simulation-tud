function [V_GO_th, V_GO_ph] = GOField(TH, PH, th_i, ph_i, F, D_r)
    theta_0 = 2 * atan(D_r / (4 * F));
    %Hyper parameters to change
    freq = 220e9;
    lambda_0 = 3e8 / freq;
    k_0 = 2 * pi / lambda_0;
    
    E_0_PW = 1;
    
    f_num = 1.8;
    D_r = 130 * lambda_0;
    
    F = D_r * f_num;
    
    R = F;
    %
    delta_n = (1 - cos(TH)) ./ (1 + cos(TH));

    amp = -2 ./ (1 + cos(TH)) .* E_0_PW;

    delta_kx = k_0 * sin(th_i) * cos(ph_i);
    delta_ky = k_0 * sin(th_i) * sin(ph_i);

    rho_fp_x = F * delta_kx / k_0;
    rho_fp_y = F * delta_ky / k_0;

    kpx = k_0 * sin(TH) .* cos(PH);
    kpy = k_0 * sin(TH) .* sin(PH);

    phase_1 = exp(-1j * (kpx * rho_fp_x + kpy * rho_fp_y));
    phase_2 = exp(-1j * ((kpx .* rho_fp_x + kpy .* rho_fp_y) .* delta_n));
    phase = phase_1 .* phase_2;

    mask = TH <= theta_0;

    E_r_GO_th = amp .* cos(PH) .* phase .* mask;
    E_r_GO_ph = -amp .* sin(PH) .* phase .* mask;

    V_GO_th = E_r_GO_th * R / exp(1j*k_0*R); 
    V_GO_ph = E_r_GO_ph * R / exp(1j*k_0*R);
end
