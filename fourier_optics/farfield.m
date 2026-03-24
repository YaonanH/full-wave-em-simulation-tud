function [Eth, Eph] = farfield(k0, R_FF, TH, PH, KZ, Gxx, Gyx, Gzx, Jx)
% Inputs:
%   k0     : free-space wavenumber
%   er_obs : relative permittivity of observation medium (use 1 for free space)
%   R_FF   : far-field radius r
%   TH,PH  : theta,phi arrays (same size)
%   Gxx,Gyx,Gzx,Jx : arrays evaluated at (kxs,kys) corresponding to TH/PH (same size)
%
% Outputs:
%   Eth,Eph : spherical components (same size)

    if ~isequal(size(TH), size(PH), size(Gxx), size(Gyx), size(Gzx), size(Jx))
        error('TH, PH, Gxx, Gyx, Gzx, and Jx must all be the same size.');
    end

    k = k0;
    common = exp(-1i*k*R_FF) ./ (2*pi*R_FF);

    Ex = 1i .* KZ .* Gxx .* Jx .* common;
    Ey = 1i .* KZ .* Gyx .* Jx .* common;
    Ez = 1i .* KZ .* Gzx .* Jx .* common;

    cth = cos(TH);  sth = sin(TH);
    cph = cos(PH);  sph = sin(PH);

    Eth = Ex .* cth .* cph + Ey .* cth .* sph - Ez .* sth;
    Eph = -Ex .* sph + Ey .* cph;
end