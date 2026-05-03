clc; clear; close all;

dx = 20e-3;
dy = 20e-3;
w = 1e-3;
l = 15e-3;

f = 10e9;
c0 = 3e8;
er = 1;
r_obs = 1;

theta_obs_deg = linspace(-90, 90, 2001);
theta_obs = deg2rad(theta_obs_deg);
phi_obs_deg = 0;
phi_obs = deg2rad(phi_obs_deg);

theta_inc_cases_deg = [0, 30];
theta_inc_cases = deg2rad(theta_inc_cases_deg);
phi_inc_deg = 0;
phi_inc = deg2rad(phi_inc_deg);

% Sweep to get AEP
theta_inc_sweep_deg = linspace(-90, 90, 2001);
theta_inc_sweep = deg2rad(theta_inc_sweep_deg);

V_tm = 1;

Mx = 10;
My = 10;
mx = -Mx:Mx;
my = -My:My;
[MX, MY] = meshgrid(mx, my);

k0 = 2 * pi * f / c0;

far_field_mag = zeros(numel(theta_inc_cases), numel(theta_obs));

for idx_case = 1:numel(theta_inc_cases)
    theta_inc = theta_inc_cases(idx_case);

    kx_inc = k0 * sin(theta_inc) * cos(phi_inc);
    ky_inc = k0 * sin(theta_inc) * sin(phi_inc);

    % Active current coefficient from the incidence direction.
    kxm = kx_inc - 2 * pi * MX / dx;
    kym = ky_inc - 2 * pi * MY / dy;

    Gxx_fs = freeSpaceGxx(k0, kxm, kym);
    basis_inc = basisFunction(k0, kxm, kym, w, l);

    v_term = cos(theta_inc) * cos(phi_inc) * V_tm / sqrt(dx * dy);
    v = v_term * basisFunction(k0, -kx_inc, -ky_inc, w, l);

    Z_terms = -(1 / (dx * dy)) * Gxx_fs .* abs(basis_inc).^2;
    Zin_active = sum(Z_terms, 'all');
    i_BF = v ./ Zin_active;

    for idx_theta = 1:numel(theta_obs)
        % Observation direction (theta, phi).
        kx_obs = k0 * sin(theta_obs(idx_theta)) * cos(phi_obs);
        ky_obs = k0 * sin(theta_obs(idx_theta)) * sin(phi_obs);

        % Windowed current spectrum:
        % Jw(kx_obs,ky_obs,kx_inc,ky_inc) =
        % iBF(kx_inc,ky_inc) * I(kx_obs) * Jt(ky_obs) * AFx * AFy.
        basis_obs = basisFunction(k0, kx_obs, ky_obs, w, l);
        AFx = arrayFactor1D(kx_obs, kx_inc, dx, Mx);
        AFy = arrayFactor1D(ky_obs, ky_inc, dy, My);
        Jw = i_BF * basis_obs .* AFx .* AFy;

        % Far field pattern:
        % e(r,theta,phi) ~ j*kz_obs*G^{ej}(kx_obs,ky_obs)*Jw*xhat*e^{-jk0r}/(2*pi*r)
        E_ff = scaledGejXColumn(er, k0, kx_obs, ky_obs) * Jw ...
            * exp(-1i * k0 * r_obs) / (2 * pi * r_obs);

        far_field_mag(idx_case, idx_theta) = norm(E_ff);
    end
end

% Active element pattern:
% sweep the incidence angle itself, as in the AEP definition.
aep_mag = zeros(size(theta_inc_sweep));

for idx_inc = 1:numel(theta_inc_sweep)
    theta_inc = theta_inc_sweep(idx_inc);

    kx_inc = k0 * sin(theta_inc) * cos(phi_inc);
    ky_inc = k0 * sin(theta_inc) * sin(phi_inc);

    kxm = kx_inc - 2 * pi * MX / dx;
    kym = ky_inc - 2 * pi * MY / dy;

    Gxx_fs = freeSpaceGxx(k0, kxm, kym);
    basis_inc = basisFunction(k0, kxm, kym, w, l);

    v_term = cos(theta_inc) * cos(phi_inc) * V_tm / sqrt(dx * dy);
    v = v_term * basisFunction(k0, -kx_inc, -ky_inc, w, l);

    Z_terms = -(1 / (dx * dy)) * Gxx_fs .* abs(basis_inc).^2;
    Zin_active = sum(Z_terms, 'all');
    i_BF = v ./ Zin_active;

    J_aep = i_BF * basisFunction(k0, kx_inc, ky_inc, w, l);
    E_aep = scaledGejXColumn(er, k0, kx_inc, ky_inc) * J_aep ...
        * exp(-1i * k0 * r_obs) / (2 * pi * r_obs);

    aep_mag(idx_inc) = norm(E_aep);
end

far_field_norm = far_field_mag / max(far_field_mag, [], 'all');
far_field_db = 20 * log10(max(far_field_norm, 1e-6));
aep_norm = aep_mag / max(aep_mag);
aep_db = 20 * log10(max(aep_norm, 1e-6));

figure('Color', 'w');
hold on; box on; grid on;
plot(theta_obs_deg, max(far_field_db(1, :), -40), 'b-', 'LineWidth', 1.3);
plot(theta_obs_deg, max(far_field_db(2, :), -40), 'r--', 'LineWidth', 1.3);
plot(theta_inc_sweep_deg, max(aep_db, -40), 'k:', 'LineWidth', 1.6);
xline(0, 'k:', 'LineWidth', 0.8);
yline(0, 'k:', 'LineWidth', 0.8);
xlim([-90, 90]);
ylim([-40, 0]);
xlabel('\theta [deg]');
ylabel('Pattern [dB]');
title(sprintf(['Far Field and AEP at 10 GHz, \\phi_{obs} = %d^\\circ, ', ...
    '\\phi_{inc} = %d^\\circ'], phi_obs_deg, phi_inc_deg));
legend('Far field, \theta_{inc} = 0^\circ', 'Far field, \theta_{inc} = 30^\circ', ...
    'AEP sweep', 'Location', 'southwest');

function Gxx = freeSpaceGxx(k0, kx, ky)
    eta0 = 120 * pi;
    kz = spectralKz(k0, kx, ky);
    Gxx = -(eta0 ./ (2 * k0 .* kz)) .* (k0^2 - kx.^2);
end

function col = scaledGejXColumn(er, k0, kx, ky)
    k = k0 * sqrt(er);
    zeta = (120 * pi) / sqrt(er);
    kz = spectralKz(k, kx, ky);

    % This is j*kz_obs*G^{ej}(kx_obs,ky_obs)*xhat, written in a form
    % that stays finite when kz_obs -> 0 at grazing observation angles.
    col = (-1j * zeta / (2 * k)) * [ ...
        k^2 - kx^2; ...
        -kx * ky; ...
        -kx * kz];
end

function af = arrayFactor1D(k_obs, k_inc, d, N)
    n = 0:(N - 1);
    af = 0;

    for idx = 1:numel(n)
        af = af + exp(1i * (k_obs - k_inc) * n(idx) * d);
    end
end

function kz = spectralKz(k0, kx, ky)
    kt2 = kx.^2 + ky.^2;
    kz = zeros(size(kt2));

    propagating = kt2 <= k0^2;
    kz(propagating) = sqrt(k0^2 - kt2(propagating));
    kz(~propagating) = -1i * sqrt(kt2(~propagating) - k0^2);
end

function basis_fun = basisFunction(k0, kx, ky, w, l)
    numerator = 2 * k0 .* (cos(kx * l / 2) - cos(k0 * l / 2));
    denominator = (k0^2 - kx.^2) .* sin(k0 * l / 2);

    Ix = numerator ./ denominator;

    singular_mask = abs(abs(kx) - k0) < 1e-10 * max(k0, 1);
    Ix(singular_mask) = l / 2;

    Jt = sinc0(ky * w / 2);
    basis_fun = Ix .* Jt;
end

function y = sinc0(x)
    y = ones(size(x));
    nz = abs(x) > 1e-12;
    y(nz) = sin(x(nz)) ./ x(nz);
end
