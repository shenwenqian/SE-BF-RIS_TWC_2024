function [g_est] = g_est_LS(H_RB, h_RU, A_B_e, A_R_e, rho, P_n, PNR_dB, T)

    [N_R, N_B] = size(H_RB);
    L = size(A_B_e,2);
    D = zeros(L,T);
    y = zeros(T,1);

    P_e = 10 ^ (PNR_dB / 10) * P_n;

    for t = 1:T
        psi = exp(1i * pi * (2 * rand(N_R,1) - ones(N_R,1)));
        f_e = randn(N_B,1) + 1i * randn(N_B,1);
        f_e = f_e / norm(f_e,'fro');
        n = (randn(1,1) + 1i * randn(1,1)) / sqrt(2) * sqrt(P_n);
        y(t) = sqrt(rho) * sqrt(P_e) * h_RU' * diag(psi) * H_RB * f_e + n;
        D(:,t) = sqrt(rho * P_e * N_B * N_R^2 / L) * ((f_e.' * conj(A_B_e)) .* (psi.' * A_R_e)).';
    end
    
    G = (conj(D) * D.') ^ (-1) * conj(D);
    g_est = G * y;
   
end

