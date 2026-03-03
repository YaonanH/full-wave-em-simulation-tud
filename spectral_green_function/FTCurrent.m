function Jx = FTCurrent(k0, er, KX, KY, l, w)

    if ~isequal(size(KX), size(KY))
        error('KX and KY must be the same size mesh grids (use meshgrid).');
    end

    keq = k0 * sqrt(er);

    x = KY .* (w/2);
    T = sinc(x/pi);

    denom = (keq^2 - KX.^2) .* sin(keq*l/2);
    L = (2*keq) .* (cos(KX*l/2) - cos(keq*l/2)) ./ denom;

    Jx = L .* T;
end