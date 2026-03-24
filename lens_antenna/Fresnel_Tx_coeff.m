function [Tau_prp, Tau_par, th_t] = Fresnel_Tx_coeff(th_i, er)
%FRESNEL_TX_COEFF Fresnel transmission from dielectric (er) to air.
    sin_th_t = sqrt(er) .* sin(th_i);
    th_t = asin(sin_th_t);

    cos_th_i = cos(th_i);
    cos_th_t = sqrt(1 - sin_th_t.^2);
    z0 = 120 * pi;
    zd = z0 / sqrt(er);

    Tau_prp = 2 * z0 * cos_th_i ./ (z0 * cos_th_i + zd * cos_th_t);
    Tau_par = 2 * z0 * cos_th_i ./ (z0 * cos_th_t + zd * cos_th_i);
end
