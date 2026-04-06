clc; clear; close all;

dx = 15e-3;
dy = 15e-3;
w = 1e-3;
l = 14e-3;

f = 10e9;
c0 = 3e8;

theta_deg = linspace(0, 90, 1000);
theta = deg2rad(theta_deg);
phi_cases_deg = [0, 90];
phi_cases = deg2rad(phi_cases_deg);
V_tm = 1;

Mx = 10;
My = 10;
mx = -Mx:Mx;
my = -My:My;
[MX, MY] = meshgrid(mx, my);

k0 = 2 * pi * f / c0;

reflection_coeff = zeros(numel(phi_cases), numel(theta));
transmission_coeff = zeros(numel(phi_cases), numel(theta));

for idx_phi = 1:numel(phi_cases)
    phi = phi_cases(idx_phi);

    for idx_theta = 1:numel(theta)
        kx_scan = k0 * sin(theta(idx_theta)) * cos(phi);
        ky_scan = k0 * sin(theta(idx_theta)) * sin(phi);

        kxm = kx_scan - 2 * pi * MX / dx;
        kym = ky_scan - 2 * pi * MY / dy;

        Gxx_fs = freeSpaceGxx(k0, kxm, kym);
        basis_fun = basisFunction(k0, kxm, kym, w, l);

        % Keep the same excitation choice already used in this file.
        v_term = V_tm / sqrt(dx * dy);
        v = v_term * basisFunction(k0, -kx_scan, -ky_scan, w, l);

        Z_terms = -(1 / (dx * dy)) * Gxx_fs .* abs(basis_fun).^2;
        Zin_active = sum(Z_terms, 'all');

        i_BF = v / Zin_active;

        Gxx_fs_scat = freeSpaceGxx(k0, kx_scan, ky_scan);
        basis_fun_scat = basisFunction(k0, kx_scan, ky_scan, w, l);
        E_scat = Gxx_fs_scat * basis_fun_scat * i_BF / (dx * dy);

        E_inc = V_tm / sqrt(dx * dy);
        Sigma_scat = E_scat / E_inc;
        Tau_scat = 1 + Sigma_scat;

        reflection_coeff(idx_phi, idx_theta) = abs(Sigma_scat);
        transmission_coeff(idx_phi, idx_theta) = abs(Tau_scat);
    end
end

for idx_phi = 1:numel(phi_cases_deg)
    figure('Color', 'w');
    hold on; box on; grid on;
    plot(theta_deg, reflection_coeff(idx_phi, :), 'k-', 'LineWidth', 1.2);
    plot(theta_deg, transmission_coeff(idx_phi, :), 'r--', 'LineWidth', 0.9);
    xline(0, 'k:', 'LineWidth', 0.8);
    xline(90, 'k:', 'LineWidth', 0.8);
    yline(0, 'k:', 'LineWidth', 0.8);
    xlim([0, 90]);
    xlabel('\theta [deg]');
    ylabel('| \Sigma |, | \Tau |');
    title(sprintf('\\phi = %d^\\circ', phi_cases_deg(idx_phi)));
    legend('Reflection', 'Transmission', 'Location', 'best');
end

function Gxx = freeSpaceGxx(k0, kx, ky)
    eta0 = 120 * pi;
    kz = spectralKz(k0, kx, ky);
    Gxx = -(eta0 ./ (2 * k0 .* kz)) .* (k0^2 - kx.^2);
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
