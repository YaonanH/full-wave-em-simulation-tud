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

z = a * sqrt(1 - (rho.^2 ./ b^2)) + c;
th_lens = atan(rho ./ z);

th_i_lens = acos((1 - eccentricity * cos(th_lens)) ./ ...
    sqrt(1 - 2 * eccentricity * cos(th_lens) + eccentricity^2));

[tau_prp_lens, tau_par_lens, th_t_lens] = Fresnel_Tx_coeff(th_i_lens, er);
tx_i_prp_ratio_lens = abs(tau_prp_lens).^2 .* (eta_d * cos(th_t_lens)) ./ (cos(th_i_lens) * eta_0);
tx_i_par_ratio_lens = abs(tau_par_lens).^2 .* (eta_d * cos(th_t_lens)) ./ (cos(th_i_lens) * eta_0);

figure;hold on; grid on;
plot(rad2deg(th_lens(:)), tx_i_prp_ratio_lens(:), 'r-', 'LineWidth', 1.2);
plot(rad2deg(th_lens(:)), tx_i_par_ratio_lens(:), 'b--', 'LineWidth', 1.2);
xlabel('\theta (deg)');
ylabel('Transmission Coefficient');
title('Fresnel Transmission Coefficients for the Lens');
legend('Parallel Polarization - TM', 'Perpendicular Polarization - TE', 'Location', 'south');
xlim([0 90]);
ylim([0 1]);

%%  Truncated lens 
D = 10 * lambda_0; % lens diameter
theta_0_truncated = deg2rad(50); % angular domain of the truncated lens
e = 1 / sqrt(er);

rho = linspace(eps, D/2 + eps, 101);
phi_lens = linspace(eps, 2*pi + eps, 101);
[RHO, PHI_lens] = meshgrid(rho, phi_lens);

r_min_truncated = (D / 2) / sin(theta_0_truncated);
a_truncated = r_min_truncated * (1 - e * cos(theta_0_truncated)) / (1 - e^2);
c_truncated = a_truncated * e;
b_truncated = sqrt(a_truncated^2 - c_truncated^2);

Z_truncated = a_truncated * sqrt(1 - (RHO.^2 ./ b_truncated^2)) + c_truncated;
TH_truncated = atan(RHO ./ Z_truncated);
r_truncated = a_truncated * (1 - e^2) ./ (1 - e * cos(TH_truncated));

[Jsx_truncated, Jsy_truncated] = LensAperture(TH_truncated, PHI_lens, r_truncated, er, u0, v0, D, theta_0_truncated);


U_truncated = sin(TH_truncated) .* cos(PHI_lens);
V_truncated = sin(TH_truncated) .* sin(PHI_lens);

Jsx_ref = max(abs(Jsx_truncated(:)));
Jsy_ref = max(abs(Jsy_truncated(:)));

Jsx_truncated_dB = 20 * log10(abs(Jsx_truncated) / Jsx_ref);
Jsy_truncated_dB = 20 * log10(abs(Jsy_truncated) / Jsy_ref);

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
xlim([0 50]);
ylim([-10 2]);

figure;
surf(U_truncated, V_truncated, Jsx_truncated_dB, 'LineStyle', 'none');
view(2); axis equal;
xlabel('U');
ylabel('V');
title('|J_{sx}| (dB)');
xlim([-1 1]);
ylim([-1 1]);
clim([-30 -20]);
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
clim([-10 0]);
colormap(jet);
colorbar;

%% Far field of the truncated lens antenna
k0 = 2 * pi / lambda_0;

phi_ff = phi_lens;

TH_ff = TH_truncated;
PH_ff = PHI_lens;
R_ff = 1;

K_xf = k0 * sin(TH_ff) .* cos(PH_ff);
K_yf = k0 * sin(TH_ff) .* sin(PH_ff);
K_zf = k0 * cos(TH_ff);

J_sx_ff = zeros(size(TH_ff));
J_sy_ff = zeros(size(TH_ff));

for m = 1:size(TH_ff, 1)
    for n = 1:size(TH_ff, 2)
        phase = exp(1j * (K_xf(m, n) * RHO .* cos(PHI_lens) + ...
            K_yf(m, n) * RHO .* sin(PHI_lens)));

        J_sx_ff(m, n) = trapz(phi_lens, trapz(rho, Jsx_truncated .* phase .* RHO, 2));
        J_sy_ff(m, n) = trapz(phi_lens, trapz(rho, Jsy_truncated .* phase .* RHO, 2));
    end
end

G_ff = EJ_SGF(1, k0, K_xf, K_yf);
common_ff = exp(-1j * k0 * R_ff) ./ (2 * pi * R_ff);

Ex_ff = 1j * K_zf .* (G_ff(:,:,1,1) .* J_sx_ff + G_ff(:,:,1,2) .* J_sy_ff) .* common_ff;
Ey_ff = 1j * K_zf .* (G_ff(:,:,2,1) .* J_sx_ff + G_ff(:,:,2,2) .* J_sy_ff) .* common_ff;
Ez_ff = 1j * K_zf .* (G_ff(:,:,3,1) .* J_sx_ff + G_ff(:,:,3,2) .* J_sy_ff) .* common_ff;

Eth_ff = Ex_ff .* cos(TH_ff) .* cos(PH_ff) + Ey_ff .* cos(TH_ff) .* sin(PH_ff) - Ez_ff .* sin(TH_ff);
Eph_ff = -Ex_ff .* sin(PH_ff) + Ey_ff .* cos(PH_ff);

Eabs_ff = sqrt(abs(Eth_ff).^2 + abs(Eph_ff).^2);
Eabs_ff_dB_norm = 20 * log10(Eabs_ff / (max(Eabs_ff(:)) + eps));
Eabs_ff_dB_norm = max(Eabs_ff_dB_norm, -40);

phi_ff_deg = rad2deg(phi_ff);
[~, idx_phi_0_lens] = min(abs(phi_ff_deg - 0));
[~, idx_phi_90_lens] = min(abs(phi_ff_deg - 90));

Eabs_ff_dB_norm_cut_phi_0 = Eabs_ff_dB_norm(idx_phi_0_lens, :);
Eabs_ff_dB_norm_cut_phi_90 = Eabs_ff_dB_norm(idx_phi_90_lens, :);

figure;
plot(rad2deg(TH_truncated(idx_phi_0_lens, :)), Eabs_ff_dB_norm_cut_phi_0, 'LineWidth', 1.5);
hold on; grid on;
plot(rad2deg(TH_truncated(idx_phi_90_lens, :)), Eabs_ff_dB_norm_cut_phi_90, 'LineWidth', 1.5);
xlabel('\theta (deg)');
ylabel('|E_{ff}| (dB)');
title('Normalized Far Field Cuts');
legend('\phi = 0 deg. cut', '\phi = 90 deg. cut', 'Location', 'northeast');

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
