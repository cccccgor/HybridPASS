function rate  = sum_rate_cal(H,F,sigma2,K)
         rate = 0;
         for k = 1 : K
             P_S = abs(H(:,k)' * F(:,k))^2;         % Signal Power
             P_I = sum(abs(H(:,k)' * F).^2) - P_S;  % Noiser Power
             rate = rate + log2(1 + P_S/(sigma2 + P_I));          % Sum-rate
         end
end
