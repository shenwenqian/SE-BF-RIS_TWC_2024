function [g_s_est] = g_est_PA(H_RB, h_RU, A_B_e, A_R_e, L_set, rho, P_n, PNR_dB, Bit, C_Ls)
    
    [N_R, N_B] = size(H_RB);
    L = size(A_B_e,2);
    L_s = length(L_set);
    y = zeros(L_s,1);

    P_e = 10 ^ (PNR_dB / 10) * P_n;

    for t = 1:L_s
        f_e = A_B_e(:,L_set(t));
        psi = N_R * conj(A_R_e(:,L_set(t)));
        n = (randn(1,1) + 1i * randn(1,1)) / sqrt(2) * sqrt(P_n);
        y(t) = sqrt(rho) * sqrt(P_e) * h_RU' * diag(psi) * H_RB * f_e + n;
    end 

    g_s_est = sqrt(L / rho / P_e / N_B / N_R^2) * y;

    if Bit ~= 0

        miu = zeros(2 ^ Bit,1);
        g_s_est_norm = g_s_est / norm(g_s_est,'fro');
        for i = 1:2 ^ Bit
            miu(i) = abs(g_s_est_norm' * C_Ls(:,i)) ^ 2;
        end
        [~, i_index] = min(miu);
        g_qua = norm(g_s_est,'fro') * C_Ls(:,i_index(1)); 
        g_s_est = g_qua;

    end

end


