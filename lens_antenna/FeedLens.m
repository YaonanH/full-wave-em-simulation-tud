function [Eth_feed, Eph_feed, Prad_feed] = FeedLens(f, er, r, TH, PH, u0, v0)
%FEEDLENS Far-field of the Gaussian lens feed inside the dielectric.
%   [Eth_feed, Eph_feed, Prad_feed] = FeedLens(f, er, r, TH, PH, u0, v0)
%   evaluates the theta/phi components shown in the lecture slide and the
%   radiated power inside the lens material.

eps_0 = 8.854187e-12;
mu_0 = 4*pi*1e-7;
zeta_0 = sqrt(mu_0/eps_0);

c0 = 3e8;
k0 = 2*pi*f/c0;
k_d = sqrt(er) * k0;
zeta_d = zeta_0 / sqrt(er);

u = sin(TH) .* cos(PH);
v = sin(TH) .* sin(PH);

amp = exp(-((u ./ u0).^2 + (v ./ v0).^2));
spherical_term = exp(-1j * k_d * r) ./ r;

Eth_feed = amp .* sin(PH) .* spherical_term;
Eph_feed = amp .* cos(PH) .* spherical_term;

Eabs_sq = abs(Eth_feed).^2 + abs(Eph_feed).^2;
U_feed = Eabs_sq .* (r.^2) / (2 * zeta_d);

theta_vec = unique(TH(1, :), "stable");
phi_vec = unique(PH(:, 1), "stable");

if isequal(size(TH), [numel(phi_vec), numel(theta_vec)]) && ...
        isequal(size(PH), [numel(phi_vec), numel(theta_vec)])
    integrand_theta = U_feed .* sin(TH);
    Prad_feed = trapz(phi_vec, trapz(theta_vec, integrand_theta, 2));
else
    warning("FeedLens:PowerIntegrationSkipped", ...
        "Prad_feed set to NaN because TH/PH are not a tensor-product mesh.");
    Prad_feed = NaN;
end
end
