clc; clear; close all;

% Assumption for the scan planes of an x-directed dipole array:
% E-plane: phi = 0 deg, D-plane: phi = 45 deg, H-plane: phi = 90 deg.

c0 = 3e8;
f0 = 10e9;
k0 = 2 * pi * f0 / c0;

theta_deg = linspace(0, 85, 681);
theta = deg2rad(theta_deg);

plane_names = {'E-plane', 'D-plane', 'H-plane'};
phi_list = [0, pi / 4, pi / 2];
plane_colors = {'k', [0.85 0.2 0.2], [0.1 0.3 0.85]};

cases(1) = struct('dx', 15e-3, 'dy', 15e-3, 'w', 1e-3, 'l', 14e-3);
cases(2) = struct('dx', 20e-3, 'dy', 20e-3, 'w', 1e-3, 'l', 14e-3);

Mx = 40;
My = 40;
mx = -Mx:Mx;
my = -My:My;
[MX, MY] = meshgrid(mx, my);

for idx_case = 1:numel(cases)
    Zin_planes = zeros(numel(phi_list), numel(theta));

    for idx_plane = 1:numel(phi_list)
        phi = phi_list(idx_plane);

        for idx_theta = 1:numel(theta)
            kx_scan = k0 * sin(theta(idx_theta)) * cos(phi);
            ky_scan = k0 * sin(theta(idx_theta)) * sin(phi);

            kxm = kx_scan + 2 * pi * MX / cases(idx_case).dx;
            kym = ky_scan + 2 * pi * MY / cases(idx_case).dy;

            Gxx_fs = freeSpaceGxx(k0, kxm, kym);
            Ix = dipoleSpectrum(kxm, k0, cases(idx_case).l);
            Jt = sinc0(kym * cases(idx_case).w / 2);

            Z_terms = -(1 / (cases(idx_case).dx * cases(idx_case).dy)) .* ...
                Gxx_fs .* abs(Ix).^2 .* abs(Jt).^2;

            Zin_planes(idx_plane, idx_theta) = sum(Z_terms, 'all');
        end
    end

    makeAngleFigure(theta_deg, Zin_planes, plane_names, plane_colors, cases(idx_case), f0);
    reportObservations(theta_deg, Zin_planes, plane_names, cases(idx_case));
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

function makeAngleFigure(theta_deg, Zin_planes, plane_names, plane_colors, params, f0)
    figure('Color', 'w');
    tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    for idx_plane = 1:numel(plane_names)
        nexttile;
        hold on; box on; grid on;
        plot(theta_deg, real(Zin_planes(idx_plane, :)), '-', ...
            'Color', plane_colors{idx_plane}, 'LineWidth', 1.3);
        plot(theta_deg, imag(Zin_planes(idx_plane, :)), '--', ...
            'Color', plane_colors{idx_plane}, 'LineWidth', 1.0);
        yline(0, 'k:', 'LineWidth', 0.8);
        xlabel('\theta [deg]');
        ylabel('Z_{in} [Ohm]');
        title(plane_names{idx_plane});
        legend('Real', 'Imaginary', 'Location', 'best');
        styleAngleAxes(gca, theta_deg);
    end

    sgtitle(sprintf(['Active impedance at %.1f GHz   |   dx = %.0f mm, dy = %.0f mm, ', ...
        'w = %.0f mm, l = %.0f mm'], ...
        f0 / 1e9, params.dx * 1e3, params.dy * 1e3, params.w * 1e3, params.l * 1e3), ...
        'FontWeight', 'bold');
end

function styleAngleAxes(ax, theta_deg)
    xlim(ax, [theta_deg(1), theta_deg(end)]);
    xticks(ax, [0 15 30 45 60 75 85]);
    ax.FontName = 'Times New Roman';
    ax.FontSize = 10;
    ax.GridColor = [0.75, 0.75, 0.75];
    ax.GridAlpha = 0.55;
    ax.LineWidth = 0.8;
end

function reportObservations(theta_deg, Zin_planes, plane_names, params)
    fprintf('\nCase: dx = %.0f mm, dy = %.0f mm, w = %.0f mm, l = %.0f mm\n', ...
        params.dx * 1e3, params.dy * 1e3, params.w * 1e3, params.l * 1e3);

    for idx_plane = 1:numel(plane_names)
        re_start = real(Zin_planes(idx_plane, 1));
        re_end = real(Zin_planes(idx_plane, end));
        im_start = imag(Zin_planes(idx_plane, 1));
        im_end = imag(Zin_planes(idx_plane, end));

        fprintf('%s: Re{Zin} %.2f -> %.2f Ohm, Im{Zin} %.2f -> %.2f Ohm for theta %.0f -> %.0f deg\n', ...
            plane_names{idx_plane}, re_start, re_end, im_start, im_end, ...
            theta_deg(1), theta_deg(end));
    end
end
