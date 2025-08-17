%% Hybrid Active and Passive Pinching Antennas
% This this the simulation codes for 
% Pinching Antenna Systems (PASS): Active, Passive, or Hybrid? 
% If you have any question with the codes, 
% please contact Kangjian Chen via kjchen@seu.edu.cn
clear all;
close all;
tic
%% PASS Setting
c = 3e8;            % Speed of light
fc = 30e9;          % Carrier Frequency
lambda = c/fc;      % Carrier Wavelength
Nt = 1000;          % Total Number of PAs
d1 = lambda/2 * 8;  % Antenna Spacing for Passive PAs
d2 = lambda/2;      % Antenna Spacing for Active PAs

spacing = lambda/2; % Minimum Spacing
Na = 24;            % Number of Active PAs
Np = 8;             % Number of Passive PAs

%% Basic Setting
K = 4;              % Number of Users
Nrf = 4;            % Number of Waveguides/RF chains
Ns = Nt/Nrf;        % Number of Antennas in each waveguide
D = (Ns-Na) * d1 + Na * d2;         % Array Aperture
L = 6;              % Number of Channel Paths
sigma2 = 1;         % Normalized Noise
%% Coordinate
% Coordinate of waveguide
Height = 3;         % Height of Waveguides        (in meters)
Width = 10;         % width of the service region (in meters)
Length = 30;        % length of the service region  (in meters)
BsAxisX = kron(linspace(Length/2/Nrf,Length-Length/2/Nrf,Nrf)',ones(Ns,1));                     % Candidate Antennas x-axis 
[BsAxisY1,x_left, x_middle, x_right] = generate_array_positions(Ns, Na, d1, d2); 
BsAxisY = kron(ones(Nrf,1),BsAxisY1'); % Candidate Antennas y-axis  (Hybrid Active Passive)
BsAxisZ = ones(Nt,1) * Height;             % Candidate Antennas z-axis 
BsAxisY2 = kron(ones(Nrf,1),(-(Ns-1)/2 : 1 : (Ns-1)/2)' * d1); % Candidate Antennas y-axis  (Only Passive)

% Communication Region
Rmin = [0,-Width/2];              % minimum communication distance [X-Axis,Y-Axis]
Rmax = [Length,Width/2];          % maximum communication distance [X-Axis,Y-Axis]


% Feasible Regions of Passive PAs
UserYRegionLeft = zeros(Nrf,1);   % Feasible Regions of Passive PAs for only Passive
UserYRegionRight = zeros(Nrf,1);  
UserYRegionLeft2 = zeros(Nrf,1);  % Active PAs divide the waveguide into two parts, their regions are denoted by UserYRegionLeft2/UserYRegionRight2 and UserYRegionLeft3/UserYRegionRight3
UserYRegionRight2 = zeros(Nrf,1);
UserYRegionLeft3 = zeros(Nrf,1);
UserYRegionRight3 = zeros(Nrf,1);
InWaveguidePositionNorm = zeros(Nrf,1);
for nn = 1 : Nrf
    UserYRegionLeft(nn) = -D/2;
    UserYRegionRight(nn) = D/2;
    UserYRegionLeft2(nn) = x_left(1);
    UserYRegionRight2(nn) = x_left(end);
    UserYRegionLeft3(nn) = x_right(1);
    UserYRegionRight3(nn) = x_right(end);
end



%% Simulation
% Channel Generation
[H,H2,ChannelPara] = ChannelGenerationContinuous(Nt,K,L,Rmin,Rmax,BsAxisX,BsAxisY,BsAxisY2,BsAxisZ,lambda);

% Noise Power
NoisePower = 10^(-11);             % linear
% Power Normarlization to prevent numerical inaccuracy due to very small values
TransmitPowerdBM = 10;
TransmitPowerLinear = 10^((TransmitPowerdBM - 30)/10);
Factor = 1/(mean(abs(H(:)).^2));
H = H * sqrt(Factor);   % Normalize the average channel power (squared magnitude) to 1
H2 = H2 * sqrt(Factor);
TransmitPowerLinear = 1/NoisePower/Factor * TransmitPowerLinear; % Adjust transmit power so that the noise power is 1

%% Hybrid Active and Passive PAs with Discrete-Position Passive PAs
[Frf1,Fbb1,Rate_Store1] =  PT_JADB(Nt,Na,Np,Nrf,sigma2,TransmitPowerLinear,K,H);
%% Hybrid Active and Passive PAs with Continuous-Position Passive PAs
Init2 = Initial_AM_JPOB(Frf1,Fbb1,Nrf,Nt,Na,Np,BsAxisY,H);      
[Frf2,Fbb2,Rate_Store2,Position2,TradeOffStore2] = AM_JPOB(Na,Nrf,sigma2,TransmitPowerLinear,K,Np,spacing,lambda,ChannelPara,BsAxisX,UserYRegionLeft2,UserYRegionRight2,UserYRegionLeft3,UserYRegionRight3,InWaveguidePositionNorm,Factor,Init2);
%% Only Active PAs 
[Frf3,Fbb3,Rate_Store3] =  PT_JADB(Nt,Na,0,Nrf,sigma2,TransmitPowerLinear,K,H);

%% Only Passive PAs with Discrete-Position Passive PAs
[Frf4,Fbb4,Rate_Store4] =  PT_JADB(Nt,0,Np,Nrf,sigma2,TransmitPowerLinear,K,H2);
%% Only Passive PAs with Continuous-Position Passive PAs
Init5 = Initial_AM_JPOB_withoutActive(Frf4,Fbb4,Nt,Nrf,Np,BsAxisY2);      
[Frf5,Fbb5,Rate_Store5,Position5,TradeOffStore5] = AM_JPOB_withoutActive(Nrf,sigma2,TransmitPowerLinear,K,Np,spacing,lambda,ChannelPara,BsAxisX,UserYRegionLeft,UserYRegionRight,InWaveguidePositionNorm,Factor,Init5);

%% Plot
plot(Rate_Store3(2:end),'-','color',[0 0.5 0],'linewidth',1.5); hold on;
plot(Rate_Store4(2:end),'k-','linewidth',1.5);
plot(Rate_Store5(2:end),'m-','linewidth',1.5);
plot(Rate_Store1(2:end),'b-','linewidth',1.5);
plot(Rate_Store2(2:end),'r-','linewidth',1.5);
xlabel('Number of iterations');
ylabel('Sum-Rate (bps/Hz)')
legend('Active PASS','Discrete-Position Passive PASS','Continuous-Position Passive PASS','Hybrid PASS with Discrete-Position Passive PAs','Hybrid PASS with Continuous-Position Passive PAs')
toc