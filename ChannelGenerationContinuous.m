function [H,H2,ChannelPara] = ChannelGenerationContinuous(Nt,K,L,Rmin,Rmax,BsAxisX,BsAxisY,BsAxisY2,BsAxisZ,lambda)
H = zeros(Nt,K);
H2 = zeros(Nt,K);

neff = 1.44;
lambdag = lambda/neff;

Rminx = Rmin(1);
Rminy = Rmin(2);
Rmaxx = Rmax(1);
Rmaxy = Rmax(2);
% 
UserAxisXStore = zeros(K,L);
UserAxisYStore = zeros(K,L);
AmplitudeStore = zeros(K,L);
for kk = 1 : K
    for ll = 1 : L
        UserAxisX = rand * (Rmaxx - Rminx) + Rminx;
        UserAxisY = rand *  (Rmaxy - Rminy) + Rminy;
        UserAxisZ = 0;
        
        FreeSpacePropogation = sqrt((UserAxisX - BsAxisX).^2 + (UserAxisY - BsAxisY).^2 + (UserAxisZ - BsAxisZ).^2);
        FreeSpacePropogation2 = sqrt((UserAxisX - BsAxisX).^2 + (UserAxisY - BsAxisY2).^2 + (UserAxisZ - BsAxisZ).^2);
        
        InWaveGuidePropogation = BsAxisY;
        InWaveGuidePropogation2 = BsAxisY2;
        
        PhasePropogation = exp(- 1j * 2 * pi * FreeSpacePropogation/lambda - 1j * 2 * pi * InWaveGuidePropogation/lambdag);
        PhasePropogation2 = exp(- 1j * 2 * pi * FreeSpacePropogation2/lambda - 1j * 2 * pi * InWaveGuidePropogation2/lambdag);

        if ll == 1 % LoS path
            Amplitude = lambda./FreeSpacePropogation/4/pi;
            Amplitude2 = lambda./FreeSpacePropogation2/4/pi;

            Factor = 1;
        else % NLoS paths
            Amplitude = lambda./FreeSpacePropogation/4/pi;
            Amplitude2 = lambda./FreeSpacePropogation/4/pi;

            Factor = (randn + 1j * randn)/sqrt(2) * 0.1;
            Amplitude = Amplitude * Factor;
            Amplitude2 = Amplitude2 * Factor;

        end
        H(:,kk) = H(:,kk) + Amplitude.* PhasePropogation;
        H2(:,kk) = H2(:,kk) + Amplitude2.* PhasePropogation2;

        % - - - - Store - - - -
        UserAxisXStore(kk,ll) = UserAxisX;
        UserAxisYStore(kk,ll) = UserAxisY;
        AmplitudeStore(kk,ll) = Factor;
    end
end

% - - - - Store - - - -
ChannelPara.UserAxisXStore = UserAxisXStore;
ChannelPara.UserAxisYStore = UserAxisYStore;
ChannelPara.AmplitudeStore = AmplitudeStore;
ChannelPara.Height = BsAxisZ;
ChannelPara.L = L;
end