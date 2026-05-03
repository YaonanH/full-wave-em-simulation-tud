clc; clear; close all;

f = 100e9;
er = 2;

r = 1;
E0 = 1;
ff_chunk_size = 128;

lambda_0 = 3e8 / f;
eta_0 = 120 * pi;
eta_d = eta_0 / sqrt(er);

a_wg = lambda_0 / 1.5;
b_wg = a_wg / 2;

theta = linspace(eps, pi/2, 181);
phi = linspace(eps, 2*pi, 181);
[TH, PH] = meshgrid(theta, phi);

[Eth_feed, Eph_feed, Prad_feed] = FeedLens_RectWG(f, er, r, TH, PH, a_wg, b_wg, E0);

U = sin(TH) .* cos(PH);
V = sin(TH) .* sin(PH);
Eabs = sqrt(abs(Eth_feed).^2 + abs(Eph_feed).^2);
Eabs_dB = 20 * log10(abs(Eabs));
Eabs_dB_norm = Eabs_dB - max(Eabs_dB(:));
Eabs_dB_norm = max(Eabs_dB_norm, -40);

phi_deg = rad2deg(phi);
theta_deg = rad2deg(theta);

% locate the indexes for phi = 0, 90, 180, 270 degrees
[~, idx_phi_0] = min(abs(phi_deg - 0));
[~, idx_phi_90] = min(abs(phi_deg - 90));
[~, idx_phi_180] = min(abs(phi_deg - 180));
[~, idx_phi_270] = min(abs(phi_deg - 270));

cut_phi_0 = Eabs_dB_norm(idx_phi_0, :);
cut_phi_90 = Eabs_dB_norm(idx_phi_90, :);
cut_phi_180 = Eabs_dB_norm(idx_phi_180, :);
cut_phi_270 = Eabs_dB_norm(idx_phi_270, :);

theta_plot = [-fliplr(theta_deg), theta_deg];
cut_phi0_phi180 = [fliplr(cut_phi_180), cut_phi_0];
cut_phi90_phi270 = [fliplr(cut_phi_270), cut_phi_90];

cut_phi0_phi180 = max(cut_phi0_phi180, -40);
cut_phi90_phi270 = max(cut_phi90_phi270, -40);

figure;
plot(theta_plot, cut_phi0_phi180, 'k-', 'LineWidth', 1.2);
hold on; grid on;
plot(theta_plot, cut_phi90_phi270, 'b--', 'LineWidth', 1.2);
xlabel('\theta (deg)');
ylabel('|E / E_{max}| (dB)');
title('Normalized farfield of the feed');
legend('\phi = 0/180 deg', '\phi = 90/270 deg', 'Location', 'south');
xlim([-90 90]);
ylim([-40 0]);

figure;
surface(U, V, Eabs_dB_norm, 'LineStyle', 'none');
axis equal;
axis([-1 1 -1 1]);
xlabel('U = sin\theta cos\phi');
ylabel('V = sin\theta sin\phi');
title('E_{feed}');
colormap(jet);
colorbar;
view(2);

%% Plot transmitted power with respect to the incident power
theta_i = linspace(eps, deg2rad(80), 181);
theta_i_deg = rad2deg(theta_i);

[tau_prp, tau_par, th_t] = Fresnel_Tx_coeff(theta_i, er);

figure; hold on; grid on;

tx_i_prp_ratio = abs(tau_prp).^2 .* (eta_d * cos(th_t)) ./ (cos(theta_i) * eta_0);
tx_i_par_ratio = abs(tau_par).^2 .* (eta_d * cos(th_t)) ./ (cos(theta_i) * eta_0);

plot(theta_i_deg, real(tx_i_prp_ratio), 'r-', 'LineWidth', 1.2);
plot(theta_i_deg, real(tx_i_par_ratio), 'b--', 'LineWidth', 1.2);
xlabel('\theta_i (deg)');
ylabel('P_t / P_{inc}');
title('Transmitted power ratio');
legend('Parallel Polarization - TM', 'Perpendicular Polarization - TE', 'Location', 'south');
xlim([0 35]);
ylim([0 1]);

%% Fresnel transmission coefficients for an untruncated elliptical fused silica lens
D_l = 18 * lambda_0; % lens diameter
eccentricity = 1 / sqrt(er);
theta_c = asin(eccentricity);
theta_0 = pi/2 - theta_c; % angular domain of the lens

r_min = (D_l / 2) / sin(theta_0);
a = r_min * (1 - eccentricity * cos(theta_0)) / (1 - eccentricity^2);
c = a * eccentricity;
b = sqrt(a^2 - c^2);

rho = linspace(eps, D_l/2, 181);

z = a * sqrt(1 - (rho.^2 ./ b^2)) + c;
th_lens = atan(rho ./ z);

th_i_lens = acos((1 - eccentricity * cos(th_lens)) ./ ...
    sqrt(1 - 2 * eccentricity * cos(th_lens) + eccentricity^2));

[tau_prp_lens, tau_par_lens, th_t_lens] = Fresnel_Tx_coeff(th_i_lens, er);
tx_i_prp_ratio_lens = abs(tau_prp_lens).^2 .* (eta_d * cos(th_t_lens)) ./ (cos(th_i_lens) * eta_0);
tx_i_par_ratio_lens = abs(tau_par_lens).^2 .* (eta_d * cos(th_t_lens)) ./ (cos(th_i_lens) * eta_0);

figure; hold on; grid on;
plot(rad2deg(th_lens(:)), tx_i_prp_ratio_lens(:), 'r-', 'LineWidth', 1.2);
plot(rad2deg(th_lens(:)), tx_i_par_ratio_lens(:), 'b--', 'LineWidth', 1.2);
xlabel('\theta (deg)');
ylabel('P_t / P_{inc}');
title('Transmitted power ratio for the Lens');
legend('Parallel Polarization - TM', 'Perpendicular Polarization - TE', 'Location', 'south');
xlim([0 70]);
ylim([0 1]);

%% Plot the equivalent aperture current distribution of an elliptical lens antenna
D = 18 * lambda_0; % lens diameter
theta_0_truncated = deg2rad(46); % angular domain of the truncated lens
e = 1 / sqrt(er);

rho = linspace(eps, D/2, 181);
phi_lens = linspace(eps, 2*pi, 181);
[RHO, PHI_lens] = meshgrid(rho, phi_lens);

r_min_truncated = (D / 2) / sin(theta_0_truncated);
a_truncated = r_min_truncated * (1 - e * cos(theta_0_truncated)) / (1 - e^2);
c_truncated = a_truncated * e;
b_truncated = sqrt(a_truncated^2 - c_truncated^2);

Z_truncated = a_truncated * sqrt(1 - (RHO.^2 ./ b_truncated^2)) + c_truncated;
TH_truncated = atan(RHO ./ Z_truncated);
r_truncated = a_truncated * (1 - e^2) ./ (1 - e * cos(TH_truncated));

[Jsx_truncated, Jsy_truncated] = LensAperture_RectWG(TH_truncated, PHI_lens, r_truncated, er, f, a_wg, b_wg, E0);

U_truncated = sin(TH_truncated) .* cos(PHI_lens);
V_truncated = sin(TH_truncated) .* sin(PHI_lens);

Jsx_truncated_dB = 20 * log10(abs(Jsx_truncated));
Jsy_truncated_dB = 20 * log10(abs(Jsy_truncated));

[~, idx_phi_0_lens] = min(abs(rad2deg(phi_lens) - 0));
[~, idx_phi_90_lens] = min(abs(rad2deg(phi_lens) - 90));

Jsy_cut_phi_0 = Jsy_truncated_dB(idx_phi_0_lens, :);
Jsy_cut_phi_90 = Jsy_truncated_dB(idx_phi_90_lens, :);

figure;
plot(rad2deg(TH_truncated(idx_phi_0_lens, :)), Jsy_cut_phi_0, 'LineWidth', 1.5);
hold on; grid on;
plot(rad2deg(TH_truncated(idx_phi_90_lens, :)), Jsy_cut_phi_90, 'LineWidth', 1.5);
xlabel('\theta (deg)');
ylabel('|J_{sy}| (dB)');
title('Not normalized cuts');
legend('\phi = 0 deg. cut', '\phi = 90 deg. cut', 'Location', 'northeast');

figure;
surf(U_truncated, V_truncated, Jsx_truncated_dB, 'LineStyle', 'none');
view(2); axis equal;
xlabel('U');
ylabel('V');
title('|J_{sx}| (dB)');
xlim([-1 1]);
ylim([-1 1]);
colormap(jet);
colorbar;

figure;
surf(U_truncated, V_truncated, Jsy_truncated_dB, 'LineStyle', 'none');
view(2); axis equal;
xlabel('U');
ylabel('V');
title('|J_{sy}| (dB)');
xlim([-1 1]);
ylim([-1 1]);
colormap(jet);
colorbar;

%% Far-field pattern
k0 = 2 * pi / lambda_0;

theta_ff = linspace(eps, pi/2 - eps, 181);
phi_ff = linspace(eps, 2*pi, 181);
phi_ff(end) = [];

[TH_ff, PH_ff] = meshgrid(theta_ff, phi_ff);
R_ff = 1;
dth_ff = theta_ff(2) - theta_ff(1);
dph_ff = phi_ff(2) - phi_ff(1);

K_xf = k0 * sin(TH_ff) .* cos(PH_ff);
K_yf = k0 * sin(TH_ff) .* sin(PH_ff);
K_zf = k0 * cos(TH_ff);

tic;
[J_sx_ff, J_sy_ff] = polar_aperture_transform_chunked( ...
    Jsx_truncated, Jsy_truncated, RHO, PHI_lens, rho, phi_lens, K_xf, K_yf, ff_chunk_size);
ff_transform_time = toc;
fprintf('Chunked far-field aperture transform time: %.2f s\n', ff_transform_time);

G_ff = EJ_SGF(1, k0, K_xf, K_yf);
[Eth_ff_x, Eph_ff_x] = farfield(k0, R_ff, TH_ff, PH_ff, K_zf, ...
    G_ff(:,:,1,1), G_ff(:,:,2,1), G_ff(:,:,3,1), J_sx_ff);
[Eth_ff_y, Eph_ff_y] = farfield(k0, R_ff, TH_ff, PH_ff, K_zf, ...
    G_ff(:,:,1,2), G_ff(:,:,2,2), G_ff(:,:,3,2), J_sy_ff);

Eth_ff = Eth_ff_x + Eth_ff_y;
Eph_ff = Eph_ff_x + Eph_ff_y;

Eabs_ff = sqrt(abs(Eth_ff).^2 + abs(Eph_ff).^2);
Eabs_ff_dB_norm = 20 * log10(Eabs_ff / (max(Eabs_ff(:)) + eps));
Eabs_ff_dB_norm = max(Eabs_ff_dB_norm, -40);

phi_ff_deg = rad2deg(phi_ff);
[~, idx_phi_0_lens] = min(abs(phi_ff_deg - 0));
[~, idx_phi_90_lens] = min(abs(phi_ff_deg - 90));

Eabs_ff_dB_norm_cut_phi_0 = Eabs_ff_dB_norm(idx_phi_0_lens, :);
Eabs_ff_dB_norm_cut_phi_90 = Eabs_ff_dB_norm(idx_phi_90_lens, :);

a_airy = D / 2;
x_airy = k0 * a_airy * sin(theta_ff);
airy = 2 * pi * a_airy^2 .* besselj(1, x_airy) ./ x_airy;
airy(x_airy == 0) = pi * a_airy^2;
airy_dB_norm = 20 * log10(abs(airy) / max(abs(airy)));
airy_dB_norm = max(airy_dB_norm, -40);

figure;
hold on; grid on;
plot(rad2deg(theta_ff), Eabs_ff_dB_norm_cut_phi_0, 'LineWidth', 1.5);
plot(rad2deg(theta_ff), Eabs_ff_dB_norm_cut_phi_90, 'LineWidth', 1.5);
plot(rad2deg(theta_ff), airy_dB_norm, 'k--', 'LineWidth', 1.2);

xlabel('\theta (deg)');
ylabel('|E_{ff}| (dB)');
title('Normalized Far Field Cuts');
legend('\phi = 0 deg. cut', '\phi = 90 deg. cut', 'Airy pattern', 'Location', 'northeast');
xlim([0 90]);
ylim([-40 0]);

U_ff = sin(TH_ff) .* cos(PH_ff);
V_ff = sin(TH_ff) .* sin(PH_ff);

figure;
surf(U_ff, V_ff, Eabs_ff_dB_norm, 'LineStyle', 'none');
view(2); axis equal;
xlabel('U');
ylabel('V');
title('|E_{ff}| (dB)');
xlim([-1 1]);
ylim([-1 1]);
clim([-40 0]);
colormap(jet);
colorbar;

%% Estimate directivity and gain of the lens antenna
[D_lens, P_rad_lens] = Directivity(Eabs_ff, TH_ff, dth_ff, dph_ff, 1, R_ff);
eta_rad_lens = P_rad_lens / Prad_feed;
G_lens = D_lens * eta_rad_lens;

D_max = max(D_lens(:));
G_max = max(G_lens(:));

fprintf('Lens radiated power: %.6g W\n', P_rad_lens);
fprintf('Radiation efficiency: %.2f%%\n', 100 * eta_rad_lens);
fprintf('Maximum directivity: %.2f dBi\n', 10 * log10(D_max));
fprintf('Maximum gain: %.2f dBi\n', 10 * log10(G_max));

function [Jx_ff, Jy_ff] = polar_aperture_transform_chunked(Jx, Jy, RHO, PHI, rho, phi, KX, KY, chunk_size)
    x_ap = RHO .* cos(PHI);
    y_ap = RHO .* sin(PHI);

    w_rho = trapz_weights(rho);
    w_phi = trapz_weights(phi);
    quad_weights = w_phi(:) * w_rho(:).';

    weighted_Jx = Jx .* RHO .* quad_weights;
    weighted_Jy = Jy .* RHO .* quad_weights;

    x_ap = x_ap(:);
    y_ap = y_ap(:);
    weighted_Jx = weighted_Jx(:);
    weighted_Jy = weighted_Jy(:);

    kx_vec = KX(:);
    ky_vec = KY(:);
    n_obs = numel(kx_vec);

    Jx_ff_vec = complex(zeros(n_obs, 1));
    Jy_ff_vec = complex(zeros(n_obs, 1));

    % Chunk the spectral points to avoid the original O(N_obs) loop over trapz calls.
    for start_idx = 1:chunk_size:n_obs
        stop_idx = min(start_idx + chunk_size - 1, n_obs);
        idx = start_idx:stop_idx;

        phase = exp(1j * (x_ap * transpose(kx_vec(idx)) + y_ap * transpose(ky_vec(idx))));
        Jx_ff_vec(idx) = phase.' * weighted_Jx;
        Jy_ff_vec(idx) = phase.' * weighted_Jy;
    end

    Jx_ff = reshape(Jx_ff_vec, size(KX));
    Jy_ff = reshape(Jy_ff_vec, size(KY));
end

function w = trapz_weights(x)
    w = zeros(size(x));

    if numel(x) == 1
        w(:) = 1;
        return;
    end

    dx = diff(x);
    w(1) = dx(1) / 2;
    w(end) = dx(end) / 2;

    if numel(x) > 2
        w(2:end-1) = (dx(1:end-1) + dx(2:end)) / 2;
    end
end
