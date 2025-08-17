function Init = Initial_AM_JPOB(Frf,Fbb,Nrf,Nt,Na,Np,BsAxisY,H)

Ns = Nt/Nrf;            % number of antennas for each waveguide
Nc = Ns - Na;
HalfN = Nc/2;           % half of candidate passive PAs


Init.Fbb = Fbb;
InitPosition = [];
InitFrf = cell(Nrf);
for kk = 1 : Nrf
    InitFrfkk = [];
    for nn = 1 : Ns
        if (nn >= HalfN + 1) && (nn <= HalfN + Na)
           InitFrfkk = [InitFrfkk;Frf((kk - 1) * Ns + nn,kk)];
        end
    end
    for nn = 1 : Ns
        if (nn <= HalfN) || (nn >= HalfN + Na + 1)
            if abs(Frf((kk - 1) * Ns + nn,kk)) > 0 
                InitPosition = [InitPosition,BsAxisY(nn)];
                InitFrfkk = [InitFrfkk;Frf((kk - 1) * Ns + nn,kk)];
            end
        end
    end
    InitFrf{kk,1} = InitFrfkk;
end

InitFrf2 = zeros((Na + Np) * Nrf,Nrf);
for kk = 1 : Nrf
 InitFrf2((kk - 1) * (Na + Np) + 1 : kk * (Na + Np),kk) =  InitFrf{kk,1} ;
end
Init.Position = InitPosition;
Init.Frf = InitFrf2;


% In fact, in our codes, the PT_JADB takes the inital Na antennas as active PAs.
% However, the AM_JPOB takes the middle Na antennas as active PAs,
% This is a conversion.
Hn = zeros(Nt,Nrf);
for kk = 1 : Nrf
    Hn((kk-1)*Ns +1 : (kk-1)*Ns + Na,:) = H((kk-1)*Ns + HalfN + 1 : (kk-1)*Ns + HalfN + Na,:);
    Hn((kk-1)*Ns +Na + 1 : (kk-1)*Ns + Na + HalfN,:) = H((kk-1)*Ns + 1 : (kk-1)*Ns + HalfN,:);
    Hn((kk-1)*Ns +Na + HalfN + 1 : (kk-1)*Ns + Ns,:) = H((kk-1)*Ns +Na + HalfN + 1 : (kk-1)*Ns + Ns,:);
end
Init.H =Hn;
end

