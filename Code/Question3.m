global C;
global numElectrons;
global x_pos_init;
global y_pos_init;
global V_x_init;
global V_y_init;
global vx
global vy;
global x;
global y;
global magThermalSpeed;
global boxWidthScaleFactor;
global boxLengthScaleFactor;
global maxTimeStep;

C.m_o = 9.10956e-31; %kg%
C.m = 0.26*C.m_o; %kg%
C.T = 300; %K%
C.k_b = 1.38064852e-23; %m^2 kg s^-2 K^-1%

meanTime = 0.2e-12; %s%
boxWidthScaleFactor = 2e-9; %nm%
boxLengthScaleFactor = 1e-9; %nm%
numElectrons = 1000;
maxTimeStep = 1000;

initPositionAndVelocity("MB", "Fragmented");

avgVelocity_init = sum(sqrt((V_x_init.^2)+(V_y_init.^2)))/numElectrons; %average of maxwell boltzmann distribution%
dt = (1/500)*100*boxWidthScaleFactor/avgVelocity_init;

dx = zeros(1, numElectrons);
dy = zeros(1, numElectrons);
vx = zeros(1, numElectrons);
vy = zeros(1, numElectrons);
x = zeros(1, numElectrons);
y = zeros(1, numElectrons);

vx = V_x_init;
vy = V_y_init;
x = x_pos_init;
y = y_pos_init;
T = zeros(1, maxTimeStep);
t_vec = linspace(1, maxTimeStep*dt, maxTimeStep);

numElectronsToTrack = 10;
electronsToTrack = randi(numElectrons, 1, numElectronsToTrack);

x_oldvals = zeros(maxTimeStep, numElectronsToTrack);
y_oldvals = zeros(maxTimeStep, numElectronsToTrack);

Pscat = 1 - exp(-dt/meanTime); %Calculate scattering probability%

avgVelocity = zeros(1, maxTimeStep);