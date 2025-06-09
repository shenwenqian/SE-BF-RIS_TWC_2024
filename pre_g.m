function [Goal] = pre_g(A_B_e, A_R_e, g_est, g_exact)                          

    Iter = 100;                                                            
    eta = 10 ^ (-2);                                                       
    [N_B, L] = size(A_B_e);
    N_R = size(A_R_e,1);
    H = sqrt(N_B * N_R ^ 2 / L) * A_R_e * diag(g_est) * A_B_e';
    psi = ones(N_R,1);
    f = ones(N_B,1) / norm(ones(N_B,1),'fro');

    eta_achi = Inf;
    Goal_achi = abs(psi.' * H * f) ^ 2;
    iter = 1;                                               

    while (eta_achi > eta) && (iter <= Iter) 

        f_last = f;
        psi_last = psi;
        Goal_last = Goal_achi;

        f = H' * conj(psi) / norm(H' * conj(psi),'fro');
        psi = conj(H * f) ./ abs(H * f);
        Goal_achi = abs(psi.' * H * f) ^ 2;

        eta_achi = Goal_achi - Goal_last;
        iter = iter + 1;

    end

    if Goal_achi < Goal_last                                              
        f = f_last;
        psi = psi_last;       
    end    

    H_exact = sqrt(N_B * N_R ^ 2 / L) * A_R_e * diag(g_exact) * A_B_e';
    Goal = abs(psi.' * H_exact * f) ^ 2;

end