function [V_a_th, V_a_ph, P_rad] = FeedField(TH, PH, D_f, k, R, d_x)
    if nargin < 6
        d_x = 0;
    end

    zeta_0 = 120 * pi;
    A_f = pi * (D_f / 2)^2;

    x = k * sin(TH) * D_f / 2;
    Airy = A_f * besselj(1, x) ./ x;
    Airy(x == 0) = A_f / 2;

    k_x = k * sin(TH) .* cos(PH);
    V_a_th = Airy .* cos(TH) .* cos(PH) .* exp(1j * k_x * d_x);
    V_a_ph = -Airy .* sin(PH) .* exp(1j * k_x * d_x);

    dth = TH(1,2) - TH(1,1);
    dph = PH(2,1) - PH(1,1);
    P_rad = sum(sum((abs(V_a_th).^2 + abs(V_a_ph).^2) .* sin(TH))) * dth * dph / (2 * zeta_0);

end
