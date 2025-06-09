function [Goal] = pre_g_s(B, L_set, g_s, g)                              

    Iter = 100;                                                            
    eta = 10 ^ (-2);                                                       
    [N_B, N_R, L] = size(B);
    L_s = length(g_s);
    L_r = L - L_s;
    L_rest = setdiff((1:1:L).',L_set); 
  
    f = ones(N_B,L_s) * conj(g_s) / norm(ones(N_B,L_s) * conj(g_s),'fro');
    psi = ones(N_R,1);

    A = zeros(N_B,L);
    for l = 1:L
        A(:,l) = B(:,:,l) * conj(psi);
    end
    A_s = A(:,L_set);
    A_r = A(:,L_rest);

    Goal_last = -Inf;
    Goal_achi = norm(A_r' * f,'fro') ^ 2 + abs(g_s.' * A_s' * f) ^ 2;
    iter = 1;                                               

    while (Goal_achi - Goal_last > eta) && (iter <= Iter) 

        A_last = A;
        f_last = f;                                                        
        Goal_last = Goal_achi;

        B_s = B(:,:,L_set);
        B_r = B(:,:,L_rest);
        Gamma_r = kron(eye(L_r), f');                             
        gamma_s = kron(g_s, f);                            
        
        Inter_1_sub = zeros(L_r,N_R);                                           
        for l_r = 1:L_r
            Inter_1_sub = Inter_1_sub + conj( Gamma_r(:,(l_r - 1) * N_B + 1:l_r * N_B) * B_r(:,:,l_r) );
        end
        Inter_1 = Inter_1_sub' * Inter_1_sub;

        Inter_2_sub = zeros(1,N_R);                                            
        for l_s = 1:L_s
            Inter_2_sub = Inter_2_sub + (gamma_s((l_s - 1) * N_B + 1:l_s * N_B,1)).' * conj(B_s(:,:,l_s));
        end
        Inter_2 = Inter_2_sub' * Inter_2_sub;
        
        J_pass = Inter_1 + Inter_2;
        eta_fp = 10^(-2);
        Iter_fp = 500;
        Obj_fp_last = 0;
        Obj_fp_achi = real(psi' * J_pass * psi);
        iter_fp = 0;
        while (iter_fp < Iter_fp) && (Obj_fp_achi - Obj_fp_last > eta_fp)
            psi_last = psi;
            Obj_fp_last = Obj_fp_achi;
            psi = exp(1i * angle(J_pass * psi));
            Obj_fp_achi = real(psi' * J_pass * psi);
            iter_fp = iter_fp + 1;
        end
        if Obj_fp_achi < Obj_fp_last
            psi = psi_last;
        end

        for l = 1:L
            A(:,l) = B(:,:,l) * conj(psi);
        end
        A_s = A(:,L_set);                                          
        A_r = A(:,L_rest);

        J_act = A_r * A_r' + A_s * conj(g_s) * g_s.' * A_s';
        [Eig_Vec, Eig_val] = eig(J_act);
        [~, max_index] = max(diag(Eig_val));
        f = Eig_Vec(:,max_index) / norm(Eig_Vec(:,max_index),'fro');                            

        Goal_achi = norm(A_r' * f,'fro') ^ 2 + abs(g_s.' * A_s' * f) ^ 2;

        iter = iter + 1;

    end

    if Goal_achi < Goal_last                                              
        A = A_last;
        f = f_last;
    end

    h = A * conj(g);
    Goal = abs(h' * f) ^ 2;
    
end