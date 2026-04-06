clc; clear; close all;

dx = 15e-3;
dy = 15e-3;
w = 1e-3;
l = 14e-3;

f = linspace(6e9, 14e9, 1000);
c0 = 3e8;

theta = 0;
phi = 0;
V_tm = 1;

Mx = 10;
My = 10;
mx = -Mx:Mx;
my = -My:My;
[MX, MY] = meshgrid(mx, my);

Zin_active = zeros(1, numel(f));
v = zeros(size(Zin_active));
i_BF = zeros(size(Zin_active));

reflection_coeff = zeros(size(Zin_active));
transmission_coeff = zeros(size(Zin_active));

for idx_f = 1:numel(f)
    k0 = 2 * pi * f(idx_f) / c0;
    kx_scan = k0 * sin(theta) * cos(phi);
    ky_scan = k0 * sin(theta) * sin(phi);

    kxm = kx_scan - 2 * pi * MX / dx;
    kym = ky_scan - 2 * pi * MY / dy;

    Gxx_fs = freeSpaceGxx(k0, kxm, kym);
    basis_fun = basisFunction(k0, kxm, kym, w, l);

    v_term = cos(theta) * cos(phi) * V_tm ./ sqrt(dx .* dy);
    Z_terms = -(1 / (dx * dy)) * Gxx_fs .* abs(basis_fun).^2;
    
    v(idx_f) = v_term .* basisFunction(k0, -kx_scan, -ky_scan, w, l);
    Zin_active(idx_f) = sum(Z_terms, 'all');

    i_BF(idx_f) = v(idx_f) ./ Zin_active(idx_f);

    % calculate the transmission and reflection coefficients 
    Gxx_fs_scat = freeSpaceGxx(k0, kx_scan, ky_scan);
    basis_fun_scat = basisFunction(k0, kx_scan, ky_scan, w, l);
    E_scat = Gxx_fs_scat * basis_fun_scat * i_BF(idx_f) ./ (dx * dy);

    E_inc = V_tm ./ sqrt(dx * dy);

    Sigma_scat = E_scat / E_inc;
    Tau_scat = 1 + Sigma_scat; 

    reflection_coeff(idx_f) = abs(Sigma_scat);
    transmission_coeff(idx_f) = abs(Tau_scat);
end



figure('Color', 'w');
hold on; box on; grid on;
plot(f, real(Zin_active), 'k-', 'LineWidth', 1.2);
plot(f, imag(Zin_active), 'k--', 'LineWidth', 0.9);
xline(10e9, 'k:', 'LineWidth', 0.8);
yline(0, 'k:', 'LineWidth', 0.8);
xlabel('Frequency [Hz]');
ylabel('Z_{in,active} [Ohm]');
title('Real and Imaginary Parts');
legend('Real', 'Imaginary', 'Location', 'northwest');

figure;
hold on; box on; grid on;
plot(f, reflection_coeff, 'k-', 'LineWidth', 1.2);
plot(f, transmission_coeff, 'k--', 'LineWidth', 0.9);
plot(f, reflection_coeff.^2 + transmission_coeff.^2, 'r-', 'LineWidth', 1.2);
xline(10e9, 'k:', 'LineWidth', 0.8);
yline(0, 'k:', 'LineWidth', 0.8);
xlabel('Frequency [Hz]');
ylabel('Coefficients');
title('Real and Imaginary Parts');
legend('Reflection', 'Transmission', 'Conservation', 'Location', 'northwest');



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
    Jt = sinc0(ky * w / 2);

    basis_fun = Ix .* Jt;
end

function y = sinc0(x)
    y = ones(size(x));
    nz = abs(x) > 1e-12;
    y(nz) = sin(x(nz)) ./ x(nz);
end
