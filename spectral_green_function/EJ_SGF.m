function [ej_SGF] = EJ_SGF(er, k0, KX, KY, sgn)
%EJ_SGF  Spectral dyadic Green's function (E due to J) on a (kx,ky) mesh.
%
%   [ej_SGF] = EJ_SGF(er, k0, KX, KY)
%   [ej_SGF] = EJ_SGF(er, k0, KX, KY, sgn)
%
% Inputs:
%   er  : relative permittivity (scalar)
%   k0  : free-space wavenumber
%   KX  : kx mesh (Ny x Nx)
%   KY  : ky mesh (Ny x Nx)
%   sgn : optional (+1 or -1) for the (±kz) terms; default +1
%
% Outputs:
%   G  : Ny x Nx x 3 x 3 dyadic, where G(:,:,i,j) = G_ij on the mesh

    if nargin < 5 || isempty(sgn)
        sgn = +1;
    end

    if ~isequal(size(KX), size(KY))
        error('KX and KY must be the same size (from meshgrid).');
    end

    % Medium wavenumber and impedance (assuming mu_r = 1)
    k    = k0 * sqrt(er);
    eta0 = 120*pi;
    zeta = eta0 / sqrt(er);

    % kz definition from your slide:
    kz = -1i .* sqrt( -(k.^2 - KX.^2 - KY.^2) );

    % Prefactor
    pref = -zeta ./ (2 .* k .* kz);

    [Ny, Nx] = size(KX);
    ej_SGF = zeros(Ny, Nx, 3, 3, 'like', kz);

    % Fill entries: G(:,:,row,col) corresponds to G_row,col
    ej_SGF(:,:,1,1) = pref .* (k.^2 - KX.^2);
    ej_SGF(:,:,1,2) = pref .* (-KX.*KY);
    ej_SGF(:,:,1,3) = pref .* (-KX.*(sgn.*kz));

    ej_SGF(:,:,2,1) = pref .* (-KX.*KY);
    ej_SGF(:,:,2,2) = pref .* (k.^2 - KY.^2);
    ej_SGF(:,:,2,3) = pref .* (-KY.*(sgn.*kz));

    ej_SGF(:,:,3,1) = pref .* (-KX.*(sgn.*kz));
    ej_SGF(:,:,3,2) = pref .* (-KY.*(sgn.*kz));
    ej_SGF(:,:,3,3) = pref .* (k.^2 - kz.^2);
end