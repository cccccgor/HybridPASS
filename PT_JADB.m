function [Frf,Fbb,Rate_Store] =  PT_JADB(Nt,Na,Np,Nrf,sigma2,P,K,H)
Ns = Nt/Nrf;            % Number of Antennas for Wach Waveguide
Nc = Ns - Na;
HalfN = Nc/2;

% Initialization
Frf = zeros(Nt,K);
for k = 1 : K
    ht = H((k-1)*Ns+1:k*Ns,k);
    for nn = 1 : Ns
        if (nn >= HalfN + 1) && (nn <= HalfN + Na)
           Frf((k-1)*Ns+nn,k) = exp(1j * angle(ht(nn)));
        else
            if abs(angle(ht(nn)))<= pi/4
                Frf((k-1)*Ns+nn,k) = 1;
            end
        end
    end
end
F = Frf;
F = F/norm(F,'fro') * sqrt(P);
rate_old = sum_rate_cal(H,F,sigma2,K);
Frf = F;
Fbb = eye(Nrf);


rho = 1;
c = 0.9;
Rate_Store = rate_old;
RateOld = 0;
for out_iter = 1 : 300   % outer loop
    for in_iter = 1 : 50 % inner loop
        %-----------------------------optimize u-----------------------------
        u = zeros(1,K);
        for k = 1 : K
            u(k) = F(:,k)' * H(:,k)/(sum(abs(H(:,k)' * F).^2) + sigma2 * norm(F,'fro')^2/P);
        end
    
        %-----------------------------optimize v-----------------------------
        v = zeros(1,K);
        for k = 1 : K
            e_k = sigma2 * abs(u(k))^2 * norm(F,'fro')^2/P + abs(1 - u(k) * H(:,k)' * F(:,k))^2 + sum(abs(u.* (H(:,k)' * F)).^2) - abs(u(k)* (H(:,k)' * F(:,k))).^2;
            v(k) = 1/e_k;
        end
    
        % -----------------------------optimize F-----------------------------
        % -----------------Compute the inverse using the Woodbury identity----------------------------
        % H_eff
        H_eff = zeros(Nt, K);
        a_vec = zeros(K, 1);
        for k = 1:K
            a_vec(k) = abs(u(k))^2 * v(k);
            H_eff(:,k) = sqrt(a_vec(k)) * H(:,k);
        end

        % alpha
        alpha = 1/(2*rho) + sigma2/P * sum(a_vec);

        temp = (eye(K) + (1/alpha) * (H_eff' * H_eff));  % size KxK
        A_inv = temp \ eye(K);

        inv_R = (1/alpha) * (eye(Nt) - (1/alpha) * H_eff * A_inv * H_eff');

        
        % -----------------Compute F----------------------------
        for k = 1 : K
            f_k_tilde = Frf * Fbb(:,k);
            eta = v(k) * conj(u(k)) * H(:,k) + f_k_tilde/2/rho;
            F(:,k) = inv_R * eta;
        end
 
        % -----------------------------Optimize B-----------------------------
        F_k_line = F;
        Frf = zeros(Nt,Nrf);
        for mm = 1 : Nrf
            LossStore = inf * ones(Ns,1);
            for nn = 1 : Ns
                tt = (mm - 1) * Ns + nn;
                ft = F(tt,:).';
                fbt = Fbb(mm,:).';
                if (nn >= HalfN + 1) && (nn <= HalfN + Na)
                    alpha = fbt' * ft/(fbt' * fbt);
                    if abs(alpha) > 1
                       Frf(tt,mm) = exp(1j * angle(alpha));
                    else
                       Frf(tt,mm) = alpha;
                    end
                else
                    Loss2 = norm(ft - fbt);
                    LossStore(nn) = Loss2;
                end
            end
            [~,Index] = sort(LossStore,'ascend');
            Frf((mm - 1) * Ns + Index(1:Np),mm) = 1;
        end
 
        % -----------------------------optimize Fbb-----------------------------
        Fbb  = inv(Frf' * Frf) * Frf' * F_k_line;
        Ratenew = RateCal(H,Frf,Fbb,P,sigma2,K);
        
        % -----------------------------inner loop stop conditions-----------------------------
        if abs(Ratenew - RateOld)< 1e-2 
            break;
        end
        RateOld = Ratenew;
    end
    [Rate_Store] = [Rate_Store,RateOld];
    rho = rho * c;
    
    % -----------------------------outer loop stop conditions-----------------------------
    if out_iter > 100 && abs(Rate_Store(out_iter) - Rate_Store(out_iter-10))<0.1
        break;
    end
end
 Fbb = Fbb/norm(Frf * Fbb,'fro') * sqrt(P);
end


function Ratenew = RateCal(H,Frf,Fbb,P,sigma2,K)

    Fbb = Fbb/norm(Frf * Fbb,'fro') * sqrt(P);
    % sum-rate calculation
    Ratenew = sum_rate_cal(H,Frf * Fbb,sigma2,K);

end