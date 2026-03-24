function [P_rx, eta_ap] = RxPower(TH, PH, V_GO_th, V_GO_ph, V_a_th, V_a_ph, P_rad)
    zeta_0 = 120 * pi;

    dth = TH(1,2) - TH(1,1);
    dph = PH(2,1) - PH(1,1);
    

    V_oc = 2 / zeta_0 * sum(sum((V_a_th .* V_GO_th + V_a_ph .* V_GO_ph) .* sin(TH))) * dth * dph;
    P_rx = abs(V_oc)^2 / (16 * P_rad);

    % Assume D_r is already defined otherwise, 
    D_r = 130 * (3e8 / 220e9); % Example value, replace with actual D_r if needed
    A_r = pi * (D_r / 2)^2;

    % Assume E_0_PW is already defined
    E_0_PW = 1; % Example value, replace with actual E_0_PW if needed
    P_inc = abs(E_0_PW)^2 * A_r / (2 * zeta_0);

    eta_ap = P_rx / P_inc;
end
