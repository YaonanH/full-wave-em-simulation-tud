clc; clear; close all;

% Lens antenna - class exercise
f = 280e9;
er = 11.9;
u0 = 0.65;
v0 = 0.65;
r = 1;

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

[~, idx_phi_0] = min(abs(phi_deg - 0));
[~, idx_phi_90] = min(abs(phi_deg - 90));
[~, idx_phi_180] = min(abs(phi_deg - 180));
[~, idx_phi_270] = min(abs(phi_deg - 270));

cut_phi_0 = Eabs_dB(idx_phi_0, :) - max(Eabs_dB(:));
cut_phi_90 = Eabs_dB(idx_phi_90, :) - max(Eabs_dB(:));
cut_phi_180 = Eabs_dB(idx_phi_180, :) - max(Eabs_dB(:));
cut_phi_270 = Eabs_dB(idx_phi_270, :) - max(Eabs_dB(:));

theta_plot = [-fliplr(theta_deg), theta_deg];
cut_phi0_phi180 = [fliplr(cut_phi_180), cut_phi_0];
cut_phi90_phi270 = [-fliplr(theta_deg), theta_deg];
cut_phi90_phi270_vals = [fliplr(cut_phi_270), cut_phi_90];

cut_phi0_phi180 = max(cut_phi0_phi180, -40);
cut_phi90_phi270_vals = max(cut_phi90_phi270_vals, -40);

figure;
plot(theta_plot, cut_phi0_phi180, 'k-', 'LineWidth', 1.2);
hold on; grid on;
plot(cut_phi90_phi270, cut_phi90_phi270_vals, 'b--', 'LineWidth', 1.2);
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
