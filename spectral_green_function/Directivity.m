function [Dir, Prad] = Directivity(E_tot, Theta, dth, dph, er, r)

    if ~isequal(size(E_tot), size(Theta))
        error('E_tot and Theta must be the same size.');
    end

    eta0 = 120*pi;
    zeta = eta0 / sqrt(er);

    % Radiation intensity
    U = (abs(E_tot).^2) .* (r.^2) ./ (2*zeta);

    % Mask out invalid points (NaN/Inf)
    valid = isfinite(U) & isfinite(Theta);

    % Radiated power (ignore invalid points)
    integrand = U .* sin(Theta);
    Prad = sum(integrand(valid), 'all') * dth * dph;

    % Directivity: define only where valid; elsewhere NaN
    Dir = NaN(size(U), 'like', U);
    Dir(valid) = 4*pi .* U(valid) ./ Prad;
end