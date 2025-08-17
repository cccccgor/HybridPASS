function [Frf,Fbb,Rate_Store,Position,TradeOffStore] = AM_JPOB_withoutActive(Nrf,sigma2,P,K,Np,spacing,lambda,ChannelPara,BsAxisX,UserYRegionLeft,UserYRegionRight,InWaveguidePositionNorm,Factor,Init)
% This function is the same as AM_JPOB but without active PAs.

neff = 1.44;
lambdag = lambda/neff;


Nt =  Np * Nrf;

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
% Passive Antenna Channel Generation
Hp = sqrt(Factor) * ChannelGenerationPassive(AntennaPostion,ChannelPara,BsAxisX,lambda,lambdag,InWaveguidePositionNorm,Nrf);


% Initialize Beamformers
Frf = Init.Frf;
Fbb = Init.Fbb;
rate_old = sum_rate_cal(Hp,Frf * Fbb,sigma2,K);


% Update F
Fp = Frf * Fbb;


Rate_Store = rate_old;
RateOld = 0;
TradeOffStore = [];
for in_iter = 1 : 200
    %-----------------------------optimize u-----------------------------
    u = zeros(1,K);
    for kk = 1 : K
        u(kk) = Fp(:,kk)' * Hp(:,kk)/(sum(abs(Hp(:,kk)' * Fp).^2) + sigma2);
    end
    if in_iter > 1
        TradeOffStore = [TradeOffStore, TradeOffCalculate(Hp,Fp,u,v,K,sigma2)];
    end
    %-----------------------------optimize v-----------------------------
    v = zeros(1,K);
    for kk = 1 : K
        e_k = sigma2 * abs(u(kk))^2  + abs(1 - u(kk) * Hp(:,kk)' * Fp(:,kk))^2 + sum(abs(u.* (Hp(:,kk)' * Fp)).^2) - abs(u(kk)* (Hp(:,kk)' * Fp(:,kk))).^2;
        v(kk) = 1/e_k;
    end
    
    if in_iter > 1
        TradeOffStore = [TradeOffStore, TradeOffCalculate(Hp,Fp,u,v,K,sigma2)];
    end
    % -----------------------------optimize yk-----------------------------
    objective_function = @(x) OnlyPassivePinchingAntennaPositionOptimizationObj(x,ChannelPara,BsAxisX,lambda,lambdag,InWaveguidePositionNorm,Nrf,u,v,Fp,Factor);
    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter');
    lb  = kron(UserYRegionLeft,ones(Np,1));
    ub = kron(UserYRegionRight,ones(Np,1));
    [AntennaPostionNew, ~] = fmincon(objective_function, AntennaPostionIni', Aneq, bneq, [], [], lb,ub, [],options);
%     [AntennaPostionNew, ~] = fmincon(objective_function, AntennaPostionIni', Aneq, bneq, [], [], lb,ub);
    Hp = sqrt(Factor) * ChannelGenerationPassive(AntennaPostionNew',ChannelPara,BsAxisX,lambda,lambdag,InWaveguidePositionNorm,Nrf);
    AntennaPostionIni = AntennaPostionNew';

 
     TradeOffStore = [TradeOffStore, TradeOffCalculate(Hp,Fp,u,v,K,sigma2)];
 

    % -----------------------------Optimize Fbb-----------------------------
    Overline_Hk = zeros(K,K);
    for kk = 1 : K
        Overline_Hk(kk,:) = u(kk) * Hp(:,kk)' * Frf;
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
    Fp = Frf * Fbb;
 
    
    TradeOffStore = [TradeOffStore, TradeOffCalculate(Hp,Fp,u,v,K,sigma2)];
   
    % -----------------------------Rate Calculation-----------------------------
    Ratenew = RateCal(Hp,Frf,Fbb,P,sigma2,K);
    [Rate_Store] = [Rate_Store,Ratenew];
    RateOld = Ratenew;
    if in_iter > 20 
       if (Rate_Store(in_iter) - Rate_Store(in_iter-20))<0.05
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