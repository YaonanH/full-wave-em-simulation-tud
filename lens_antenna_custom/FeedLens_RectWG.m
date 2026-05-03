function [Eth, Eph, Prad] = FeedLens_RectWG(f, er, r, TH, PH, a, b, E0)

    c0 = 3e8;
    k0 = 2*pi*f/c0;
    k  = sqrt(er)*k0;
    eta = 120*pi/sqrt(er);

    % First-order air-to-dielectric transmission
    Gamma = (1 - sqrt(er))/(1 + sqrt(er));
    T = 1 + Gamma;
    

    % Far-zone variables
    X = (k*a/2).*sin(TH).*cos(PH);
    Y = (k*b/2).*sin(TH).*sin(PH);

    FX = cos(X)./(X.^2 - (pi/2)^2);
    FY = sin(Y)./Y;

    % Avoid numerical 0/0
    FY(abs(Y) < 1e-12) = 1;
    FX(abs(abs(X) - pi/2) < 1e-12) = -1/pi;

    C = 1j*a*b*k*(T*E0).*exp(-1j*k*r)./(2*pi*r);

    Eth = -(pi/2).*C.*sin(PH).*FX.*FY;
    Eph = -(pi/2).*C.*cos(TH).*cos(PH).*FX.*FY;

    U = (abs(Eth).^2 + abs(Eph).^2).*r.^2./(2*eta);

    dth = TH(1,2) - TH(1,1);
    dph = PH(2,1) - PH(1,1);

    Prad = sum(sum(U.*sin(TH)))*dth*dph;

end