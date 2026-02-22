function Jx = FTCurrent(k0, er, KX, KY, l, w)
%FTCurrent  Fourier transform of x-directed current density Jx(KX,KY)
%
%   Jx = FTCurrent(k0, er, KX, KY, l, w)
%
% Inputs:
%   KX, KY : kx/ky mesh grids (Ny x Nx), e.g. [KX,KY] = meshgrid(kx,ky)
%
% Based on:
%   Jx^FT(kx,ky) = L(kx) * T(ky)
%   T(ky) = sinc(ky*w/2)                      (unnormalized sinc = sin(x)/x)
%   L(kx) = 2*keq * (cos(kx*l/2) - cos(keq*l/2)) ...
%           / ((keq^2 - kx^2) * sin(keq*l/2))
%
% Here: keq = k0*sqrt(er) (mu_r assumed 1).
%
% Notes:
% - MATLAB's sinc(x) is normalized: sinc(x)=sin(pi*x)/(pi*x)
%   so sin(x)/x is sinc(x/pi).

    % Validate mesh sizes
    if ~isequal(size(KX), size(KY))
        error('KX and KY must be the same size mesh grids (use meshgrid).');
    end

    % Equivalent wavenumber in medium
    keq = k0 * sqrt(er);

    % ----- T(KY): unnormalized sinc(x)=sin(x)/x -----
    x = KY .* (w/2);
    T = sinc(x/pi);

    % ----- L(KX) -----
    denom = (keq^2 - KX.^2) .* sin(keq*l/2);

    % numerical guard around denom ~ 0 (elementwise)
    tol = 1e-12 * max(1, abs(keq^2 * sin(keq*l/2)));
    denom(abs(denom) < tol) = tol;

    L = (2*keq) .* (cos(KX*l/2) - cos(keq*l/2)) ./ denom;

    % ----- Product (same size as KX, KY) -----
    Jx = L .* T;
end