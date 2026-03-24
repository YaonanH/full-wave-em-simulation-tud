clc; clear; close all;

% Lens antenna - class exercise
f = 280e9;
er = 11.9;
u0 = 0.65;
v0 = 0.65;
r = 1;

lambda_0 = 3e8 / f;
eta_0 = 120 * pi;
eta_d = eta_0 / sqrt(er);

theta = linspace(eps, pi/2 - eps, 1000);
phi = linspace(0, 2*pi, 1000);
[TH, PH] = meshgrid(theta, phi);

[Eth_feed, Eph_feed, Prad_feed] = FeedLens(f, er, r, TH, PH, u0, v0);

U = sin(TH) .* cos(PH);
V = sin(TH) .* sin(PH);
Eabs = sqrt(abs(Eth_feed).^2 + abs(Eph_feed).^2);
Eabs_dB = 20 * log10(abs(Eabs));
Eabs_dB_norm = Eabs_dB - max(Eabs_dB(:));
Eabs_dB_norm = max(Eabs_dB_norm, -40);

phi_deg = rad2deg(phi);
theta_deg = rad2deg(theta);

%locate the indexes for phi = 0, 90, 180, 270 degrees
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
ylabel('|E| / E_{max} (dB)');
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
clim([-40 0]);
view(2);

fprintf('Radiated power inside the lens: %.6g W (for unit field scaling)\n', Prad_feed);

%% Plot transmitted power with respect to the incident power
theta_i = linspace(eps, deg2rad(20), 1000);
theta_i_deg = rad2deg(theta_i);

[tau_prp, tau_par, th_t] = Fresnel_Tx_coeff(theta_i, er);

figure;hold on; grid on;

tx_i_prp_ratio = abs(tau_prp).^2 .* (eta_d * cos(th_t)) ./ (cos(theta_i) * eta_0);
tx_i_par_ratio = abs(tau_par).^2 .* (eta_d * cos(th_t)) ./ (cos(theta_i) * eta_0);

plot(theta_i_deg, tx_i_prp_ratio, 'r-', 'LineWidth', 1.2);
plot(theta_i_deg, tx_i_par_ratio, 'b--', 'LineWidth', 1.2);
xlabel('\theta (deg)');
ylabel('Transmission Coefficient');
title('Fresnel Transmission Coefficients');
legend('Parallel Polarization - TM', 'Perpendicular Polarization - TE', 'Location', 'south');
xlim([0 20]);
ylim([0 1]);

%% Plot power ratio for a lens
D_l = 10 * lambda_0; % lens diameter
eccentricity = 1 / sqrt(er);
theta_c = asin(eccentricity);
theta_0 = pi/2 - theta_c; % angular domain of the lens

r_min = (D_l / 2) / sin(theta_0);
a = r_min * (1 - eccentricity * cos(theta_0)) / (1 - eccentricity^2);
c = a * eccentricity;
b = sqrt(a^2 - c^2);

rho = linspace(eps, D_l/2, 1000);
phi_lens = linspace(eps, 2*pi, 1000);
[RHO, PHI_lens] = meshgrid(rho, phi_lens);

z = a * sqrt(1 - (RHO.^2 ./ b^2)) + c;
TH_lens = atan(RHO ./ z);

TH_i_lens = acos((1 - eccentricity * cos(TH_lens)) ./ ...
    sqrt(1 - 2 * eccentricity * cos(TH_lens) + eccentricity^2));

[tau_prp_lens, tau_par_lens, th_t_lens] = Fresnel_Tx_coeff(TH_i_lens, er);
tx_i_prp_ratio_lens = abs(tau_prp_lens).^2 .* (eta_d * cos(th_t_lens)) ./ (cos(TH_i_lens) * eta_0);
tx_i_par_ratio_lens = abs(tau_par_lens).^2 .* (eta_d * cos(th_t_lens)) ./ (cos(TH_i_lens) * eta_0);

figure;hold on; grid on;
plot(rad2deg(TH_lens(:)), tx_i_prp_ratio_lens(:), 'r-', 'LineWidth', 1.2);
plot(rad2deg(TH_lens(:)), tx_i_par_ratio_lens(:), 'b--', 'LineWidth', 1.2);
xlabel('\theta (deg)');
ylabel('Transmission Coefficient');
title('Fresnel Transmission Coefficients for the Lens');
legend('Parallel Polarization - TM', 'Perpendicular Polarization - TE', 'Location', 'south');
xlim([0 90]);
ylim([0 1]);

