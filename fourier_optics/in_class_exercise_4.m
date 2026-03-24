%% Initialize parameters
clc; clear; close all;

freq = 300e9;
lambda_0 = 3e8 / freq;
k_0 = 2 * pi / lambda_0;

E_0_PW = 1;

f_num = 2;
D_r = 100 * lambda_0;
D_f = 4 * lambda_0;

% theta_0 = asin(1/(2*f_num));
% F = D_r ./ (2 * tan(theta_0));
F = D_r * f_num;


R = F;

th_i = deg2rad(1.13);
ph_i = 0;

%% Plot GO field and feed field on the FO sphere
theta = linspace(0, pi/2, 501);
phi = linspace(0, 2*pi, 501);

[TH, PH] = meshgrid(theta, phi);

[V_GO_th, V_GO_ph] = GOField(TH, PH, th_i, ph_i, F, D_r);

d_x = 2*lambda_0*f_num;
[V_a_th, V_a_ph, P_rad] = FeedField(TH, PH, D_f, k_0, R, d_x);

e_r_GO_th =  V_GO_th .* exp(1j * k_0 * R) ./ R;
e_r_GO_ph =  V_GO_ph .* exp(1j * k_0 * R) ./ R;

E_a_tx_th =  V_a_th .* exp(1j * k_0 * R) ./ R;
E_a_tx_ph =  V_a_ph .* exp(1j * k_0 * R) ./ R;

e_r_GO_1 = sqrt(abs(V_GO_th(1, :)).^2 + abs(V_GO_ph(1, :)).^2);
e_r_GO_2 = sqrt(abs(V_GO_ph(250, :)).^2 + abs(V_GO_th(250, :)).^2);
E_a_tx_1 = sqrt(abs(E_a_tx_th(1, :)).^2 + abs(E_a_tx_ph(1, :)).^2);
E_a_tx_2 = sqrt(abs(E_a_tx_ph(250, :)).^2 + abs(E_a_tx_th(250, :)).^2);


figure
hold on
grid on
plot(rad2deg(theta), 20 * log10(e_r_GO_1 / max(e_r_GO_1)), 'r', 'LineWidth', 1.2)
plot(rad2deg(theta), 20 * log10(e_r_GO_2 / max(e_r_GO_2)), 'k--', 'LineWidth', 1.2)
plot(rad2deg(theta), 20 * log10(E_a_tx_1 / max(E_a_tx_1)), 'k', 'LineWidth', 1.2)
plot(rad2deg(theta), 20 * log10(E_a_tx_2 / max(E_a_tx_2)), 'b--', 'LineWidth', 1.2)
xlabel('\theta (deg)')
ylabel('|E / E_{max}| (dB)')
title('Normalized GO field and Feed field on the FO sphere')
legend('|e_{GO}^{r}|, \phi = 0^\circ and \phi = 180^\circ', ...
       '|e_{GO}^{r}|, \phi = 90^\circ and \phi = 270^\circ', ...
       '|E_{a}^{tx}|, \phi = 0^\circ and \phi = 180^\circ', ...
       '|E_{a}^{tx}|, \phi = 90^\circ and \phi = 270^\circ', ...
       'Location', 'south')
xlim([-15 15])
ylim([-15 0])


[P_rx, eta_ap] = RxPower(TH, PH, V_GO_th, V_GO_ph, V_a_th, V_a_ph, P_rad);

%% Plot received power pattern
theta_rx = linspace(0, pi/2, 501);
phi_rx = linspace(0, 2 * pi, 501);
[TH_rx, PH_rx] = meshgrid(theta_rx, phi_rx);

d_x = 2*lambda_0*f_num;
[V_a_th_rx, V_a_ph_rx, P_rad_rx] = FeedField(TH_rx, PH_rx, D_f, k_0, R, d_x);

th_i_sweep = linspace(-deg2rad(3), deg2rad(3), 301);
P_rx_cut = zeros(size(th_i_sweep));

for n = 1:length(th_i_sweep)
    th_i = th_i_sweep(n);
    ph_i = 0;
    [V_GO_th_rx, V_GO_ph_rx] = GOField(TH_rx, PH_rx, th_i, ph_i, F, D_r);
    [P_rx, ~] = RxPower(TH_rx, PH_rx, V_GO_th_rx, V_GO_ph_rx, V_a_th_rx, V_a_ph_rx, P_rad_rx);
    P_rx_cut(n) = P_rx;
end

figure
plot(rad2deg(th_i_sweep), 10 * log10(P_rx_cut / max(P_rx_cut)), 'k', 'LineWidth', 1.2)
grid on
xlabel('\theta_i (deg)')
ylabel('Normalized pattern (dB)')
title('Pattern in Reception')
xlim([-3 3])
ylim([-40 0])

%% Plot power patter for displaced field
