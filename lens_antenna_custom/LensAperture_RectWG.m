function [Jx, Jy] = LensAperture_RectWG(th, ph, r, er, f, a, b, E0)
%LENSAPERTURE_RECTWG Equivalent aperture currents using rectangular waveguide TE10 feed.

    c0 = 3e8;
    k0 = 2*pi*f/c0;
    kd = sqrt(er)*k0;

    e = 1 / sqrt(er);
    eta0 = 120*pi;

    % Rectangular waveguide TE10 far-field variables
    X = (kd*a/2).*sin(th).*cos(ph);
    Y = (kd*b/2).*sin(th).*sin(ph);

    FX = cos(X)./(X.^2 - (pi/2)^2);
    FY = sin(Y)./Y;

    % Avoid numerical singularities
    FY(abs(Y) < 1e-12) = 1;
    FX(abs(abs(X) - pi/2) < 1e-12) = -1/pi;

    % First-order air-to-dielectric transmission at feed/lens interface
    T_in = 2/(1 + sqrt(er));
    

    % Full far-field amplitude constant 
    % (remove .*exp(-1j*kd*r)./r)
    C_modified = 1j*a*b*kd*(T_in*E0) ./(2*pi);

    % Rectangular waveguide TE10 feed field components
    Ei_th = -(pi/2).*C_modified.*sin(ph).*FX.*FY;
    Ei_ph = -(pi/2).*C_modified.*cos(th).*cos(ph).*FX.*FY;

    den = sqrt(1 + e^2 - 2 * e .* cos(th));
    cos_th_i = (1 - e .* cos(th)) ./ den;
    th_i = acos(cos_th_i);

    [Tau_prp, Tau_par, th_t] = Fresnel_Tx_coeff(th_i, er);

    S_f = sqrt((cos(th_t) ./ cos(th_i)) .* ((e .* cos(th) - 1) ./ (e - cos(th))));

    Jx = -(2 / eta0) .* (Tau_par .* Ei_th .* cos(ph) - Tau_prp .* Ei_ph .* sin(ph)) .* S_f ./ r;
    Jy = -(2 / eta0) .* (Tau_par .* Ei_th .* sin(ph) + Tau_prp .* Ei_ph .* cos(ph)) .* S_f ./ r;


end