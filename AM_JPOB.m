function [Frf,Fbb,Rate_Store,Position,TradeOffStore] =  AM_JPOB(Na,Nrf,sigma2,P,K,Np,spacing,lambda,ChannelPara,BsAxisX,UserYRegionLeft2,UserYRegionRight2,UserYRegionLeft3,UserYRegionRight3,InWaveguidePositionNorm,Factor,Init)
neff = 1.44;
lambdag = lambda/neff;



Nt = (Na + Np) * Nrf;
Ns = Nt/Nrf;

% inequality constraints
Aneq = zeros((Np - 1) * Nrf,Np * Nrf);
for ii = 1 : Nrf
    for jj = 1 : Np -1
        index = (ii - 1) * Np +jj;
        row = (ii - 1) * (Np -1) + jj;
        Aneq(row,index) = 1;
        Aneq(row,index + 1) = -1;
    end
end
bneq = -ones((Np - 1) * Nrf,1) * spacing;


% Position Initialization
AntennaPostion =  Init.Position;
AntennaPostionIni = AntennaPostion;
H = Init.H;
Hr = H;
% Boundary Generation
lb  = zeros(Nrf * Np,1);
ub = zeros(Nrf * Np,1);
for ii = 1 : Nrf * Np
    jj = ceil(ii/Np);
    if (AntennaPostionIni(ii) >=  UserYRegionLeft2(jj)) && (AntennaPostionIni(ii) <=  UserYRegionRight2(jj))
        lb(ii) = UserYRegionLeft2(jj);
        ub(ii) = UserYRegionRight2(jj);
    else
        lb(ii) = UserYRegionLeft3(jj);
        ub(ii) = UserYRegionRight3(jj);
    end
end



% Passive Antenna Channel Generation
Hp = sqrt(Factor) * ChannelGenerationPassive(AntennaPostion,ChannelPara,BsAxisX,lambda,lambdag,InWaveguidePositionNorm,Nrf);
% Active Antenna Channel Generation
Ha = zeros(Na * Nrf, Nrf);
Nsr = size(H,1)/Nrf;
for kk = 1 : Nrf
    Ha((kk-1)*Na +1 : kk*Na,:) = Hr((kk-1)*Nsr +1 : (kk-1)*Nsr + Na,:);
end

% Channel Update
H = zeros(Nt,Nrf);
for kk = 1 : Nrf
    H((kk-1)*Ns +1 : (kk-1)*Ns + Na,:) = Ha((kk-1)*Na +1 : kk*Na,:);
    H((kk-1)*Ns +1 + Na :  kk * Ns,:) = Hp((kk-1)*Np +1 :  kk * Np,:);
end

% Initialize Beamformers
Frf = Init.Frf;
Fbb = Init.Fbb;
rate_old = sum_rate_cal(H,Frf * Fbb,sigma2,K);


% Update F
F = Frf * Fbb;
% Update Fa 和 Fp
Fa = zeros(Na * Nrf,Nrf);
Fp = zeros(Np * Nrf,Nrf);
for kk = 1 : Nrf
    Fa((kk-1) * Na + 1 :  kk * Na,:) = F((kk-1) * Ns + 1 : (kk-1) * Ns + Na,:);
    Fp((kk-1) * Np + 1 :  kk * Np,:) = F((kk-1) * Ns + Na + 1 : kk * Ns,:);
end

Rate_Store = rate_old;
RateOld = 0;
TradeOffStore = [];  % Calculate the Objective Tradeoff
for in_iter = 1 : 200
    %-----------------------------optimize u-----------------------------
    u = zeros(1,K);
    for kk = 1 : K
        u(kk) = F(:,kk)' * H(:,kk)/(sum(abs(H(:,kk)' * F).^2) + sigma2);
    end
    
    if in_iter > 1
        TradeOffStore = [TradeOffStore, TradeOffCalculate(H,F,u,v,K,sigma2)];
    end
    %-----------------------------optimize v-----------------------------
    v = zeros(1,K);
    for kk = 1 : K
        e_k = sigma2 * abs(u(kk))^2  + abs(1 - u(kk) * H(:,kk)' * F(:,kk))^2 + sum(abs(u.* (H(:,kk)' * F)).^2) - abs(u(kk)* (H(:,kk)' * F(:,kk))).^2;
        v(kk) = 1/e_k;
    end
    
    if in_iter > 1
        TradeOffStore = [TradeOffStore, TradeOffCalculate(H,F,u,v,K,sigma2)];
    end
    % -----------------------------optimize PA Postions-----------------------------
    objective_function = @(x) PassivePinchingAntennaPositionOptimizationObj(x,ChannelPara,BsAxisX,lambda,lambdag,InWaveguidePositionNorm,Nrf,u,v,Ha,Fa,Fp,Factor);
    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');
    [AntennaPostionNew, ~] = fmincon(objective_function, AntennaPostionIni', Aneq, bneq, [], [], lb,ub,[],options);
    Hp = sqrt(Factor) * ChannelGenerationPassive(AntennaPostionNew',ChannelPara,BsAxisX,lambda,lambdag,InWaveguidePositionNorm,Nrf);
    AntennaPostionIni = AntennaPostionNew';
    % Update Channels
    H = zeros(Nt,Nrf);
    for kk = 1 : Nrf
        H((kk-1)*Ns +1 : (kk-1)*Ns + Na,:) = Ha((kk-1)*Na +1 : kk*Na,:);
        H((kk-1)*Ns +1 + Na :  kk * Ns,:) = Hp((kk-1)*Np +1 :  kk * Np,:);
    end
    
 
     TradeOffStore = [TradeOffStore, TradeOffCalculate(H,F,u,v,K,sigma2)];
 
    
    % ---------------------------Optimize Frf----------------------------
    % Obtain Frfb
    Frfb = zeros(Np * Nrf,Nrf);
    for kk = 1 : Nrf
        Frfb((kk-1) * Np + 1 : kk * Np,kk) = ones(Np,1); 
    end
    
    
    Mu = zeros(K,K);
    for kk = 1 : K
        for ii = 1 : K
            if kk == ii
                hpk = Hp(:,kk);
                fpi = Fp(:,ii);
                Mu(kk,ii) = u(kk) * hpk' * fpi - 1;
            else
                hpk = Hp(:,kk);
                fpi = Fp(:,ii);
                Mu(kk,ii) = u(kk) * hpk' * fpi;
            end
        end
    end

    Q = zeros(Na * Nrf,Na * Nrf);
    q = zeros(Na * Nrf,1);
    for kk = 1 : K
        for ii = 1 : K
            hk = Ha(:,kk);
            wi = Fbb(:,ii);
            h_ki = hk .* conj(repelem(wi, Na));
            Q = Q + v(kk) * abs(u(kk))^2 * (h_ki * h_ki');
            q = q + v(kk) * conj(u(kk)) * h_ki * Mu(kk,ii);
        end
    end
    % QCQP
     cvx_begin quiet
        variable frf(Na * Nrf,1) complex   
        expression Frfa(Na * Nrf, Nrf)     
        for kk = 1:Nrf
            row_idx = (kk-1)*Na + (1:Na);    
            Frfa(row_idx, kk) = frf((kk-1)*Na + (1:Na)); 
        end

        minimize real( quad_form(frf, Q) + 2 * real(q' * frf) )
        subject to
            for nn = 1 : Na * Nrf
                square_abs(frf(nn)) <= 1;
            end
            sum_square_abs(Frfa * Fbb) + sum_square_abs(Frfb * Fbb) <= P;
%             norm(Frfa * Fbb,'fro') <= sqrt(P - norm(Frfb * Fbb,'fro')^2);
    cvx_end
    % Reconstruct Frf
    for kk = 1 : Nrf
        for nn = 1 : Na
            Frf((kk-1) * Ns + nn,kk) = frf((kk-1) * Na + nn);
        end
    end
    % Reconstruct F;
    F = Frf * Fbb;
    % Update Fa 和 Fp
    Fa = zeros(Na * Nrf,Nrf);
    Fp = zeros(Np * Nrf,Nrf);
    for kk = 1 : Nrf
        Fa((kk-1) * Na + 1 :  kk * Na,:) = F((kk-1) * Ns + 1 : (kk-1) * Ns + Na,:);
        Fp((kk-1) * Np + 1 :  kk * Np,:) = F((kk-1) * Ns + Na + 1 : kk * Ns,:);
    end
    
    TradeOffStore = [TradeOffStore, TradeOffCalculate(H,F,u,v,K,sigma2)];

    % -----------------------------Optimize Fbb-----------------------------
    Overline_Hk = zeros(K,K);
    for kk = 1 : K
        Overline_Hk(kk,:) = u(kk) * H(:,kk)' * Frf;
    end
    T = zeros(K,K);
    S = zeros(K,K);
    for kk = 1 : K
        T = T + v(kk) * Overline_Hk(kk,:)' * Overline_Hk(kk,:);
        S(:,kk) = v(kk) * Overline_Hk(kk,:)';
    end
    overlineP = kron(eye(K),T);
    overlinep = S(:);
    overlineFrf = kron(eye(K),Frf);
    cvx_begin quiet
        variable fbb(K^2,1) complex
        minimize(( quad_form(fbb, overlineP) - 2 * real(overlinep' * fbb) ) )
        subject to
            norm(overlineFrf * fbb) <= sqrt(P);
    cvx_end
    % Reconstruct Fbb
    Fbb = reshape(fbb,K,K);
    % Reconstruct F;
    F = Frf * Fbb;
    % Update Fa 和 Fp
    Fa = zeros(Na * Nrf,Nrf);
    Fp = zeros(Np * Nrf,Nrf);
    for kk = 1 : Nrf
        Fa((kk-1) * Na + 1 :  kk * Na,:) = F((kk-1) * Ns + 1 : (kk-1) * Ns + Na,:);
        Fp((kk-1) * Np + 1 :  kk * Np,:) = F((kk-1) * Ns + Na + 1 : kk * Ns,:);
    end
    
    TradeOffStore = [TradeOffStore, TradeOffCalculate(H,F,u,v,K,sigma2)];
   
    % -----------------------------Rate Calculation-----------------------------
    Ratenew = RateCal(H,Frf,Fbb,P,sigma2,K);
    [Rate_Store] = [Rate_Store,Ratenew];
    RateOld = Ratenew;
    if in_iter > 40 
       if (Rate_Store(in_iter) - Rate_Store(in_iter-40))<0.05
            break;
       end
    end
end
Position = AntennaPostionNew;
end


function Ratenew = RateCal(H,Frf,Fbb,P,sigma2,K)

    Fbb = Fbb/norm(Frf * Fbb,'fro') * sqrt(P);
    % Rate Calculation
    Ratenew = sum_rate_cal(H,Frf * Fbb,sigma2,K);

end

function tradeoff = TradeOffCalculate(H,F,u,v,K,sigma2)
tradeoff = 0;

for kk = 1 : K
    ek = abs(u(kk))^2 * sigma2;
    for ii = 1 : K
       if  ii ~= kk
           ek = ek + abs(u(kk) * H(:,kk)' * F(:,ii))^2;
       else
           ek = ek + abs(u(kk) * H(:,kk)' * F(:,ii) - 1)^2;
       end
    end
    tradeoff = tradeoff + v(kk) * ek - log(v(kk));
end

end