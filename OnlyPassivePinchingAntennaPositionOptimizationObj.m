function   objective_function = OnlyPassivePinchingAntennaPositionOptimizationObj(x,ChannelPara,BsAxisX,lambda,lambdag,InWaveguidePositionNorm,Nrf,u,v,Fp,Factor) 


x = x';
objective_function = 0;
Hp = sqrt(Factor) * ChannelGenerationPassive(x,ChannelPara,BsAxisX,lambda,lambdag,InWaveguidePositionNorm,Nrf);

for kk = 1 : Nrf
    ek = abs(1 - u(kk) * Hp(:,kk)' * Fp(:,kk))^2 +...
        sum(abs(u.* (Hp(:,kk)' * Fp) ).^2)   - ...
        abs(u(kk)* (Hp(:,kk)' * Fp(:,kk))).^2;
    objective_function = objective_function + v(kk) * ek;
end

end
