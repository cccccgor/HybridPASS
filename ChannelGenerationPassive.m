function Hp = ChannelGenerationPassive(AntennaPostion,ChannelPara,BsAxisX,lambda,lambdag,InWaveguidePositionNorm,Nrf)
UserAxisXStore = ChannelPara.UserAxisXStore;
UserAxisYStore = ChannelPara.UserAxisYStore;
AmplitudeStore = ChannelPara.AmplitudeStore;
Height = ChannelPara.Height;

[K,L] = size(UserAxisXStore);
Npt = numel(AntennaPostion);
Hp = zeros(Npt,K);

Nsp = Npt/Nrf;
Nt = numel(BsAxisX);
Ns = Nt/Nrf;

BsAxisXN = kron(BsAxisX(1:Ns:end),ones(Nsp,1))';
for kk = 1 : K
   for ll = 1 : L
       UserAxisX = UserAxisXStore(kk,ll);
       UserAxisY = UserAxisYStore(kk,ll);
       UserAxisZ = 0;
       BsAxisY = AntennaPostion;
       BsAxisZ = Height(1);
       Factor = AmplitudeStore(kk,ll);

       FreeSpacePropogation = sqrt((UserAxisX - BsAxisXN).^2 + (UserAxisY - BsAxisY).^2 + (UserAxisZ - BsAxisZ).^2);
%        InWaveGuidePropogation =  AntennaPostion - kron(InWaveguidePositionNorm',ones(1,Nt/Nrf));
       InWaveGuidePropogation =  AntennaPostion;
       PhasePropogation = exp(- 1j * 2 * pi * FreeSpacePropogation/lambda - 1j * 2 * pi * InWaveGuidePropogation/lambdag);
       Amplitude = lambda./FreeSpacePropogation/4/pi;
       
       Hp(:,kk) = Hp(:,kk) + Factor * Amplitude' .*  PhasePropogation.';
   end
end
end