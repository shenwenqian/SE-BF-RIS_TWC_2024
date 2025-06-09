function [L_set_full] = L_sel(B, L_RB, L_RU, L_s_min)    

    [N_B, N_R, L] = size(B);
    Iter = 100;                                                           
    eta = 10 ^ (-2);  

    L_set_full = zeros(L - 1,L - 1);                                            
    L_set = (1:1:L).';

    for L_s = L:-1:L_s_min + 1                                                  
        
        L_r = L - L_s;
        L_rest = setdiff((1:1:L).',L_set);
        
        V = ones(N_B,L_s) / norm(ones(N_B,L_s),'fro');      
        psi = ones(N_R,1);

        A = zeros(N_B,L);
        for l = 1:1:L       
            A(:,l) = B(:,:,l) * conj(psi);                             
        end
        A_s = A(:,L_set);
        A_r = A(:,L_rest);

        Q = (L + L_RB + L_RU - 3) / (L - 1);

        Goal_last = -Inf;
        if L_s ~= L
            Goal_achi = (4 - 2 * Q) * norm(diag(A_s' * V),'fro') ^ 2 + Q * abs(trace(A_s' * V)) ^ 2 + ...
                        Q * norm(A_s' * V,'fro') ^ 2 + Q * norm(A_r' * V,'fro') ^ 2;
        else
            Goal_achi = (4 - 2 * Q) * norm(diag(A_s' * V),'fro') ^ 2 + Q * abs(trace(A_s' * V)) ^ 2 + ...
                        Q * norm(A_s' * V,'fro') ^ 2;
        end        

        iter = 1; 

        while (Goal_achi - Goal_last > eta) && (iter < Iter) 

            A_last = A;
            V_last = V;                                                
            Goal_last = Goal_achi;

            B_s = B(:,:,L_set);
            B_r = B(:,:,L_rest);
            Upsilon_s = kron(eye(L_s),V');  
            Upsilon_r = kron(eye(L_r),V');  
        
            Inter_1 = zeros(N_R,N_R);                                            
            for l_s = 1:L_s
                Inter_1_sub = V(:,l_s).' * conj(B_s(:,:,l_s));
                Inter_1 = Inter_1 + Inter_1_sub' * Inter_1_sub;
            end
        
            Inter_2_sub = zeros(1,N_R);
            for l_s = 1:L_s
                Inter_2_sub = Inter_2_sub + V(:,l_s).' * conj(B_s(:,:,l_s));
            end
            Inter_2 = Inter_2_sub' * Inter_2_sub;

            Inter_3_sub = zeros(L_s * L_s,N_R);
            for l_s = 1:L_s
                Inter_3_sub = Inter_3_sub + conj( Upsilon_s(:,(l_s - 1) * N_B + 1:l_s * N_B) ) * conj( B_s(:,:,l_s) );
            end
            Inter_3 = Inter_3_sub' * Inter_3_sub;

            Inter_4_sub = zeros(L_s * L_r,N_R);
            for l_r = 1:L_r
                Inter_4_sub = Inter_4_sub + conj( Upsilon_r(:,(l_r - 1) * N_B + 1:l_r * N_B) ) * conj( B_r(:,:,l_r) );
            end
            Inter_4 = Inter_4_sub' * Inter_4_sub;
        
            if L_s ~= L
                J_pass = (4 - 2 * Q) * Inter_1 + Q * Inter_2 + Q * Inter_3 + Q * Inter_4;
            else
                J_pass = (4 - 2 * Q) * Inter_1 + Q * Inter_2 + Q * Inter_3;
            end
                                 
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

            for l = 1:1:L       
                A(:,l) = B(:,:,l) * conj(psi);                             
            end
            A_s = A(:,L_set);
            A_r = A(:,L_rest);     

            Inter_5 = [];
            for l_s = 1:L_s
                Inter_5 = blkdiag(Inter_5, A_s(:,l_s) * A_s(:,l_s)');
            end
        
            if L_s ~= L
                J_act = (4 - 2 * Q) * Inter_5 + Q * vec(A_s) * (vec(A_s))' + ...
                        Q * (kron(eye(L_s), A_s'))' * kron(eye(L_s), A_s') + ...
                        Q * (kron(eye(L_s), A_r'))' * kron(eye(L_s), A_r');
            else 
                J_act = (4 - 2 * Q) * Inter_5 + Q * vec(A_s) * (vec(A_s))' + ...
                        Q * (kron(eye(L_s), A_s'))' * kron(eye(L_s), A_s');
            end 

            [Eig_Vec, Eig_val] = eig(J_act);
            [~, max_index] = max(diag(Eig_val));
            Max_Eig_Vec = Eig_Vec(:,max_index);                            
            V = reshape(Max_Eig_Vec,N_B,L_s) / norm(Max_Eig_Vec,'fro');

            if L_s ~= L
                Goal_achi = (4 - 2 * Q) * norm(diag(A_s' * V),'fro') ^ 2 + Q * abs(trace(A_s' * V)) ^ 2 + ...
                            Q * norm(A_s' * V,'fro') ^ 2 + Q * norm(A_r' * V,'fro') ^ 2;
            else
                Goal_achi = (4 - 2 * Q) * norm(diag(A_s' * V),'fro') ^ 2 + Q * abs(trace(A_s' * V)) ^ 2 + ...
                            Q * norm(A_s' * V,'fro') ^ 2;
            end 

            iter = iter + 1;

        end

        if Goal_achi < Goal_last       
            A = A_last;
            A_s = A(:,L_set);
            A_r = A(:,L_rest);             
            V = V_last;
        end

        Z = zeros(L_s, 1);                                         
        
        for l_s = 1:L_s

            Inter_6 = 0;
            Inter_7 = 0;
            for i = 1:L_s
                if i ~= l_s
                    Inter_6 = Inter_6 + A_s(:,i)' * V(:,i);
                    Inter_7 = Inter_7 + V(:,l_s)' * A_s(:,l_s) * A_s(:,i)' * V(:,i);
                end
            end

            if L_s ~= L                                                       
                Z(l_s) = 4 * V(:,l_s)' * A_s(:,l_s) * A_s(:,l_s)' * V(:,l_s) + ...
                         2 * Q * real( V(:,l_s)' * A_s(:,l_s) * Inter_6 ) + ...
                         2 * Q * real( Inter_7 ) + Q * norm(A_r' * V(:,l_s),'fro') ^ 2; 
            else                                                              
                Z(l_s) = 4 * V(:,l_s)' * A_s(:,l_s) * A_s(:,l_s)' * V(:,l_s) + ...
                         2 * Q * real( V(:,l_s)' * A_s(:,l_s) * Inter_6 ) + ...
                         2 * Q * real( Inter_7 );
            end


        end    

        [~, min_index] = min(Z);                                  
        L_set = setdiff(L_set,L_set(min_index(1)));
        L_set_full(1:L_s - 1,L_s - 1) = L_set;

    end

end
