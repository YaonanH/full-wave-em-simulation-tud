function [Jx, Jy] = LensAperture(th, ph, r, er, u0, v0, D, th_0)
%LENSAPERTURE Equivalent aperture currents of a truncated elliptical lens.
    e = 1 / sqrt(er);
    eta0 = 120 * pi;

    u = sin(th) .* cos(ph);
    v = sin(th) .* sin(ph);
    amp = exp(-((u ./ u0).^2 + (v ./ v0).^2));

    Ei_th = amp .* sin(ph) ./ r;
    Ei_ph = amp .* cos(ph) ./ r;

    den = sqrt(1 + e^2 - 2 * e .* cos(th));
    cos_th_i = (1 - e .* cos(th)) ./ den;
    th_i = acos(cos_th_i);

    [Tau_prp, Tau_par, th_t] = Fresnel_Tx_coeff(th_i, er);

    S_f = sqrt((cos(th_t) ./ cos(th_i)) .* ((e .* cos(th) - 1) ./ (e - cos(th))));

    Jx = -(2 / eta0) .* (Tau_par .* Ei_th .* cos(ph) - Tau_prp .* Ei_ph .* sin(ph)) .* S_f ./ r;
    Jy = -(2 / eta0) .* (Tau_par .* Ei_th .* sin(ph) + Tau_prp .* Ei_ph .* cos(ph)) .* S_f ./ r;

end
