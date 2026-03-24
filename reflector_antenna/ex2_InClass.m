clc; clear; close all;

% Reflector antenna - Matlab instruction 2
eps_0 = 8.854187e-12;
mu_0 = 4*pi*1e-7;
zeta = sqrt(mu_0/eps_0);

c = 3e8;
f = 300e9;
lambda = c/f;
Df = 5*lambda;
a = Df/2;
k = 2*pi/lambda;
Jf = 1;     % constant amplitude

theta_feed = linspace(eps, pi/2-eps, 1000);
phi_feed = linspace(eps, 2*pi-eps, 1000);
r_feed = 1;

[TH_feed, PH_feed] = meshgrid(theta_feed, phi_feed);
J_FT = 2*pi*a^2*besselj(1, k*a*sin(TH_feed))./(k*a*sin(TH_feed));
% J_FT = airy(k*a*sin(TH_feed));
Eff_th = -2j*k*zeta*Jf*J_FT.*cos(TH_feed).*sin(PH_feed)*exp(-1j*k*r_feed)/(4*pi*r_feed);
Eff_ph = -2j*k*zeta*Jf*J_FT.*cos(PH_feed)*exp(-1j*k*r_feed)/(4*pi*r_feed);

E_tot = zeros([size(Eff_ph), 3]);
E_tot(:, :, 1) = Eff_th;
E_tot(:, :, 2) = Eff_ph;

U = sin(TH_feed).*cos(PH_feed);
V = sin(TH_feed).*sin(PH_feed);

E_ff = sqrt(abs(Eff_th).^2 + abs(Eff_ph).^2);
Eff_maxdB = max(max(20*log10(abs(E_ff))));
E_ff_dB_norm  = 20*log10(abs(E_ff)) - Eff_maxdB;
figure
surface(U, V, E_ff_dB_norm , 'linestyle' , 'none' )
axis equal
colormap("jet")
colorbar
clim([-40 0])

[~, idx_phi_90] = min(abs(phi_feed - pi/2));
figure
hold on; grid on;
plot(rad2deg(theta_feed), 20*log10(abs(Eff_ph(1, :))) - Eff_maxdB, "-k", "LineWidth", 2)
plot(rad2deg(theta_feed), 20*log10(abs(Eff_th(idx_phi_90, :))) - Eff_maxdB, "--r", "LineWidth", 2)
xlabel("\theta [deg]"); ylabel("|E| normalized [dB]");
title("E^{ff}_{\phi, feed}")
legend("\phi_{feed}=0", "\phi_{feed}=\pi/2")
ylim([-40, 0])
%%
focal = 3;
D = 4;
f2 = focal/D;
theta_0 = 2*acot(4*f2);
U_feed = 1/(2*zeta)* abs(E_ff).^2 * r_feed^2;

dth = abs(TH_feed(1, 1) - TH_feed(1, 2));
dph = abs(PH_feed(1, 1) - PH_feed(2, 1));
numerator = U_feed.*sin(TH_feed);
numerator(TH_feed > theta_0) = 0;
eta_spillOver = (sum(sum(numerator))*dth*dph)/(sum(sum(U_feed.*sin(TH_feed)))*dth*dph);
fprintf("The spillover efficiency is: %.2f%% \n", eta_spillOver*100)

%%
close all; clc;
theta = linspace(0, theta_0, 1000);
[TH, PH_feed] = meshgrid(theta, phi_feed);
r_prime = 2*focal./(1+cos(TH));
rho_prime = r_prime.*sin(TH);
phi_prime = PH_feed;

drho = abs(rho_prime(1, 1) - rho_prime(1, 2));
dphi = abs(phi_prime(1, 1) - phi_prime(2, 1));

Eff_th = -2j*k*zeta*Jf*J_FT .* cos(TH).*sin(PH_feed) * exp(-1j*k*r_feed)/(4*pi*r_feed);
Eff_ph = -2j*k*zeta*Jf*J_FT .* cos(PH_feed)          * exp(-1j*k*r_feed)/(4*pi*r_feed);

Ea_rho = -Eff_th*exp(-1j*k*2*focal).*cos(TH/2).^2/focal;
Ea_phi = -Eff_ph*exp(-1j*k*2*focal).*cos(TH/2).^2/focal;

Ea_x = Ea_rho.*cos(phi_prime) - Ea_phi.*sin(phi_prime);
Ea_y = Ea_rho.*sin(phi_prime) + Ea_phi.*cos(phi_prime);
E_a = zeros([size(Ea_x), 3]);
E_a(:, :, 1) = Ea_x;
E_a(:, :, 2) = Ea_y;

n = zeros([size(TH), 3]); 
n(:,:,3) = ones(size(TH));
M = cross(E_a, n, 3);
M_x = squeeze(M(:, :, 1));
M_y = squeeze(M(:, :, 2));

M = sqrt(abs(M(:, :, 1)).^2 + abs(M(:, :, 2)).^2 + abs(M(:, :, 3)).^2);

U = rho_prime.*cos(phi_prime);
V = rho_prime.*sin(phi_prime);

figure
hold on
surface(U, V, 20*log10(abs(M_x))-max(max(20*log10(abs(M_x)))) , 'linestyle' , 'none' ) ;
title("|M_x|")
colormap("jet")
colorbar
clim([-40 0])
axis([ -2 2 -2 2])

figure
hold on
surface(U, V, 20*log10(abs(M_y)) - max(max(20*log10(abs(M_y)))) , 'linestyle' , 'none' ) ;
axis equal
axis([ -2 2 -2 2])
title("|M_y|")
colormap("jet")
colorbar
clim([-80 -40])

%%
integrand_num = sqrt(abs(E_a(:, :, 1)).^2 + abs(E_a(:, :, 2)).^2 + abs(E_a(:, :, 3)).^2);
numerator = sum(sum(integrand_num.*rho_prime));
integrand_den = sqrt(abs(E_a(:, :, 1)).^2 + abs(E_a(:, :, 2)).^2 + abs(E_a(:, :, 3)).^2);
denominator = sum(sum(abs(integrand_den).^2.*rho_prime));
A = pi*(D/2)^2;
eta_taper = 1/A * (abs(numerator*dphi*drho).^2)/(denominator*drho*dphi)

%
int_x = sum(sum(E_a(:,:,1) .* rho_prime)) * drho * dphi;
int_y = sum(sum(E_a(:,:,2) .* rho_prime)) * drho * dphi;
int_z = sum(sum(E_a(:,:,3) .* rho_prime)) * drho * dphi;

numerator_sq = abs(int_x)^2 + abs(int_y)^2 + abs(int_z)^2;

% Denominator: integral of |E|^2 * rho'
integrand_den = (abs(E_a(:,:,1)).^2 + abs(E_a(:,:,2)).^2 + abs(E_a(:,:,3)).^2);
denominator = sum(sum(integrand_den .* rho_prime)) * drho * dphi;

% Taper efficiency
A = pi * (D/2)^2;
eta_taper = (1/A) * numerator_sq / denominator

%% 
% z = 
% J_s = 1/zeta * cross(n, cross(n, E_a, 3), 3);
% E_far = zeros(size(rho_prime));
% 
% % th_ff = linspace(eps, pi/2-eps, 1000);
% ph_ff = linspace(eps, 2*pi-eps, 1000);
% th_ff = linspace(0, 5*lambda/D, 1000);
% % ph_ff = [0 pi/2];
% [TH_ff, PH_ff] = meshgrid(th_ff, ph_ff);
% 
% Kx = k*sin(TH_ff).*cos(PH_ff);
% Ky = k*sin(TH_ff).*sin(PH_ff);
% Kz = k*cos(TH_ff);
% J_ft = J_s.*exp(1j*Kx.*rho_prime.*cos(phi_prime)).*exp(1j*Ky.*rho_prime.*sin(phi_prime)).*phi_prime*drho*dphi;
% 
% GF = EJ_SGF(zeta, k, Kx, Ky);
% Gxx = squeeze(GF(1, 1, :, :));
% Gyx = squeeze(GF(1, 2, :, :));
% Gzx = squeeze(GF(1, 3, :, :));
% R_ff = 1;
% [Eff_th, Eff_ph] = FarField(k, Kz, R_ff, TH_ff, PH_ff, Gxx, Gyx, Gzx, J_ft);
% 
% E_ff = sqrt(abs(Eff_th).^2 + abs(Eff_ph).^2);
% E_ff_max = max(max(20*log10(abs(E_ff))));
% 
% [~, idx_phi_90] = min(abs(ph_ff - pi/2));
% 
% figure
% hold on; grid on;
% plot(th_ff, 20*log10(abs(E_ff(1, :))) - E_ff_max, "-k", "LineWidth", 2)
% plot(th_ff, E_ff(idx_phi_90, :) - E_ff_max, "--r", "LineWidth", 2)
% for i = 1:length(th_ff)
%     for j = 1:length(ph_ff)
% 
%     end
% end

%%
% Surface current (PMC image: factor 2)
J_s = 1/zeta * cross(n, cross(n, E_a, 3), 3);
Js_x = 2 * squeeze(J_s(:,:,1));
Js_y = 2 * squeeze(J_s(:,:,2));

% Far field grid
% ph_ff = linspace(eps, 2*pi-eps, 100);
ph_ff = [0 pi/2];
th_ff = linspace(0, 5*lambda/D, 100);


% Pre-allocate
Eth = zeros(length(ph_ff), length(th_ff));
Eph = zeros(length(ph_ff), length(th_ff));

R_ff = 1;

for ii = 1:length(th_ff)
    for jj = 1:length(ph_ff)
        
        % Current observation point
        th  = th_ff(ii);
        ph  = ph_ff(jj);
        
        kx = k*sin(th)*cos(ph);
        ky = k*sin(th)*sin(ph);
        kz = k*cos(th);
        
        % Phase term over aperture [N_ap_phi x N_ap_rho]
        phase = exp(1j*(kx*rho_prime.*cos(phi_prime) + ky*rho_prime.*sin(phi_prime)));
        
        % Integrate over aperture
        JFT_x = sum(sum(Js_x .* phase .* rho_prime)) * drho * dphi;
        JFT_y = sum(sum(Js_y .* phase .* rho_prime)) * drho * dphi;
        
        % Green's function (scalar kx, ky for this observation point)
        GF = EJ_SGF(zeta, k, kx, ky);
        Gxx = GF(1,1); Gyx = GF(2,1); Gzx = GF(3,1);
        Gxy = GF(1,2); Gyy = GF(2,2); Gzy = GF(3,2);
        
        % Far field contributions
        [Eth_x, Eph_x] = FarField(k, kz, R_ff, th, ph, Gxx, Gyx, Gzx, JFT_x);
        [Eth_y, Eph_y] = FarField(k, kz, R_ff, th, ph, Gxy, Gyy, Gzy, JFT_y);
        
        Eth(jj, ii) = Eth_x + Eth_y;
        Eph(jj, ii) = Eph_x + Eph_y;
      
    end
end

% Total field
E_ff = sqrt(abs(Eth).^2 + abs(Eph).^2);
E_ff_max = max(max(20*log10(abs(E_ff))));

[~, idx_phi_0]  = min(abs(ph_ff - 0));
[~, idx_phi_90] = min(abs(ph_ff - pi/2));

figure; hold on; grid on;
plot(rad2deg(th_ff), 20*log10(abs(E_ff(idx_phi_0,:)))  - E_ff_max, '-k', 'LineWidth', 2)
plot(rad2deg(th_ff), 20*log10(abs(E_ff(idx_phi_90,:))) - E_ff_max, '--r', 'LineWidth', 2)
xlabel('\theta [deg]'); ylabel('|E| normalized [dB]')
legend('\phi=0°', '\phi=90°'); ylim([-40 0])