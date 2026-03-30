clc; clear; close all;

dx = 15e-3;
dy = 15e-3;
w = 3e-3;
l = 14e-3;

f = linspace(5e9, 15e9, 1000);
c0 = 3e8;

theta_cases = [0, 30 * pi / 180];
phi_cases = [0, 0];

Mx = 40;
My = 40;
mx = -Mx:Mx;
my = -My:My;
[MX, MY] = meshgrid(mx, my);

idx_m0 = find(mx == 0);
idx_n0 = find(my == 0);

Zin_active = zeros(numel(theta_cases), numel(f));
Zin_fundamental = zeros(numel(theta_cases), numel(f));
Zin_higher_order = zeros(numel(theta_cases), numel(f));

for idx_case = 1:numel(theta_cases)
    theta = theta_cases(idx_case);
    phi = phi_cases(idx_case);

    for idx_f = 1:numel(f)
        k0 = 2 * pi * f(idx_f) / c0;
        kx_scan = k0 * sin(theta) * cos(phi);
        ky_scan = k0 * sin(theta) * sin(phi);

        kxm = kx_scan + 2 * pi * MX / dx;
        kym = ky_scan + 2 * pi * MY / dy;

        Gxx_fs = freeSpaceGxx(k0, kxm, kym);
        Ix = dipoleSpectrum(kxm, k0, l);
        Jt = sinc0(kym * w / 2);

        Z_terms = -(1 / (dx * dy)) * Gxx_fs .* abs(Ix).^2 .* abs(Jt).^2;

        Zin_fundamental(idx_case, idx_f) = Z_terms(idx_n0, idx_m0);
        Zin_active(idx_case, idx_f) = sum(Z_terms, 'all');
        Zin_higher_order(idx_case, idx_f) = ...
            Zin_active(idx_case, idx_f) - Zin_fundamental(idx_case, idx_f);
    end
end

for idx_case = 1:numel(theta_cases)
    makeCaseFigure( ...
        f, ...
        Zin_active(idx_case, :), ...
        Zin_higher_order(idx_case, :), ...
        Zin_fundamental(idx_case, :), ...
        theta_cases(idx_case), ...
        phi_cases(idx_case));
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

function Ix = dipoleSpectrum(kx, k0, l)
    numerator = 2 * k0 .* (cos(kx * l / 2) - cos(k0 * l / 2));
    denominator = (k0^2 - kx.^2) .* sin(k0 * l / 2);

    Ix = numerator ./ denominator;

    singular_mask = abs(abs(kx) - k0) < 1e-10 * max(k0, 1);
    Ix(singular_mask) = l / 2;
end

function y = sinc0(x)
    y = ones(size(x));
    nz = abs(x) > 1e-12;
    y(nz) = sin(x(nz)) ./ x(nz);
end

function makeCaseFigure(f, Zin_total, Zin_higher, Zin_fundamental, theta, phi)
    figure('Color', 'w');
    tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    hold on; box on; grid on;
    plot(f, real(Zin_total), 'k-', 'LineWidth', 1.2);
    plot(f, imag(Zin_total), 'k--', 'LineWidth', 0.9);
    xline(10e9, 'k:', 'LineWidth', 0.8);
    yline(0, 'k:', 'LineWidth', 0.8);
    xlabel('Frequency [Hz]');
    ylabel('Z_{in,active} [Ohm]');
    title('Real and Imaginary Parts');
    legend('Real', 'Imaginary', 'Location', 'northwest');
    styleAxes(gca, f);

    nexttile;
    hold on; box on; grid on;
    plot(f, real(Zin_higher), 'r-', 'LineWidth', 1.2);
    plot(f, imag(Zin_higher), 'r--', 'LineWidth', 0.9);
    xline(10e9, 'k:', 'LineWidth', 0.8);
    yline(0, 'k:', 'LineWidth', 0.8);
    xlabel('Frequency [Hz]');
    ylabel('Z_{higher} [Ohm]');
    title('Higher-Order Modes');
    legend('Real', 'Imaginary', 'Location', 'northwest');
    styleAxes(gca, f);

    nexttile;
    hold on; box on; grid on;
    plot(f, real(Zin_fundamental), 'b-', 'LineWidth', 1.2);
    plot(f, imag(Zin_fundamental), 'b--', 'LineWidth', 0.9);
    xline(10e9, 'k:', 'LineWidth', 0.8);
    yline(0, 'k:', 'LineWidth', 0.8);
    xlabel('Frequency [Hz]');
    ylabel('Z_{fundamental} [Ohm]');
    title('Fundamental Mode');
    legend('Real', 'Imaginary', 'Location', 'northwest');
    styleAxes(gca, f);

    sgtitle(sprintf('Theta = %.4g; Phi = %.4g;', theta, phi), ...
        'FontWeight', 'bold');
end

function styleAxes(ax, f)
    xlim(ax, [f(1), f(end)]);
    xticks(ax, [5e9, 10e9, 15e9]);
    ax.FontName = 'Times New Roman';
    ax.FontSize = 10;
    ax.GridColor = [0.75, 0.75, 0.75];
    ax.GridAlpha = 0.55;
    ax.LineWidth = 0.8;
    ax.XAxis.Exponent = 10;
end
