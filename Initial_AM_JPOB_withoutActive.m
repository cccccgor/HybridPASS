function Init = Initial_AM_JPOB_withoutActive(Frf,Fbb,Nt,Nrf,Np,BsAxisY)
Ns = Nt/Nrf;
Init.Fbb = Fbb;
InitPosition = [];
InitFrf = cell(Nrf);
for kk = 1 : Nrf
    InitFrfkk = [];
    for nn = (kk - 1) * Ns + 1 : kk * Ns
        if abs(Frf(nn,kk)) > 0 
           InitPosition = [InitPosition,BsAxisY(nn)];
           InitFrfkk = [InitFrfkk;Frf(nn,kk)];
        end
    end
    InitFrf{kk,1} = InitFrfkk;
end

InitFrf2 = zeros( Np  * Nrf,Nrf);
for kk = 1 : Nrf
 InitFrf2((kk - 1) * Np + 1 : kk * Np,kk) =  InitFrf{kk,1} ;
end
Init.Position = InitPosition;
Init.Frf = InitFrf2;
end

