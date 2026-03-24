function [Eth_feed, Eph_feed, Prad_feed] = FeedLens(f, er, r, TH, PH, u0, v0)
    k0 = 2*pi*f/3e8;
    kd = sqrt(er) * k0;
    eta_d = 120*pi / sqrt(er);

    u = sin(TH) .* cos(PH);
    v = sin(TH) .* sin(PH);
    amp = exp(-((u ./ u0).^2 + (v ./ v0).^2));
    phase = exp(-1j * kd * r) ./ r;

    Eth_feed = amp .* sin(PH) .* phase;
    Eph_feed = amp .* cos(PH) .* phase;

    U = (abs(Eth_feed).^2 + abs(Eph_feed).^2) .* r.^2 / (2 * eta_d);

    dth = TH(1,2) - TH(1,1);
    dph = PH(2,1) - PH(1,1);
    % Theta from 0 to pi, Phi from 0 to 2*pi

    Prad_feed = sum(sum((U) .* sin(TH))) * dth * dph;
end
