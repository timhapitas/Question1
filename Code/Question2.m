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

initPositionAndVelocity("MB", "Uniform");

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

numElectronsToTrack = 15;
electronsToTrack = randi(numElectrons, 1, numElectronsToTrack);

x_oldvals = zeros(maxTimeStep, numElectronsToTrack);
y_oldvals = zeros(maxTimeStep, numElectronsToTrack);

Pscat = 1 - exp(-dt/meanTime); %Calculate scattering probability%

avgVelocity = zeros(1, maxTimeStep);
figure;
for t = 1:maxTimeStep
    
    x_previous = x;
    y_previous = y;
    
    vx_previous = vx;
    vy_previous = vy;
    
    %Check boundary conditions of electrons
    
    toReflect = find((y_previous > 100*boxLengthScaleFactor) | (y_previous < 0));
    toShiftRight = find(x_previous < 0);
    toShiftLeft = find(x_previous > 100*boxWidthScaleFactor);
    
    if ~isempty(toReflect)
        vy_previous(toReflect) = -vy_previous(toReflect);
        vy = vy_previous;
    end
    
    if ~isempty(toShiftRight)
        x_previous(toShiftRight) = 100*boxWidthScaleFactor; 
    end
    
    if ~isempty(toShiftLeft)
        x_previous(toShiftLeft) = 0;
       
    end
    
    %Check if any electrons scatter. If they do, re-thermalize all that
    %have scattered (sample new speeds from maxwell boltzmann distribution)
    
    scatteringProbabilites = rand(1, numElectrons);
    electronsToScatter = find(scatteringProbabilites <= Pscat);
    [vx(electronsToScatter), vy(electronsToScatter)] = thermalize(length(electronsToScatter));
    
    dx = vx * dt;
    dy = vy * dt;
    
    x = x_previous + dx;
    y = y_previous + dy;
    
    x_oldvals(t,:) = x_previous(electronsToTrack);
    y_oldvals(t,:) = y_previous(electronsToTrack);
    
    %compute running temperature
    avgVelocity(t) = sum(sqrt((vx.^2)+(vy.^2)))/numElectrons; %average of maxwell boltzmann distribution%
    
    % plot code goes here*
    
    plot(x, y, 'b.');
    axis([0 100*boxWidthScaleFactor 0 100*boxLengthScaleFactor]);
    pause(0.0001);
    
end

%Trajectory of select electrons over entire simulation duration
figure;
title('Particle Trajectories for 15 Randomly Selected Electrons');
for pltCnt = 1:numElectronsToTrack
    
   x_diffs = diff(x_oldvals(:,pltCnt));
   x_diffs(length(x_diffs) + 1) = 0;
   
   rejections = find(abs(x_diffs) >= 100e-9);
   x_oldvals(rejections,pltCnt) = NaN;
    
   plot(x_oldvals(:,pltCnt), y_oldvals(:,pltCnt));
   hold on;
end
xlabel('Horizontal Location (m)');
ylabel('Vertical Location (m)');
axis([0 100*boxWidthScaleFactor 0 100*boxLengthScaleFactor]);
hold off;

%Temperature plot over time

T = ((avgVelocity.^2).*(C.m))./(pi.*(C.k_b)); 

figure;
plot(t_vec, T, 'r');
title("Semiconductor Temperature Vs Time");
xlabel("Time (s)");
ylabel("Temperature (K)")
grid on;

% Histogram plots of velocity components
figure;
histogram(V_x_init, 50);
title("Histogram of x-Velocities for Electrons in the Semiconductor");
xlabel("x-Velocity (m/s)");
ylabel("Number of Electrons per Bin");

figure;
histogram(V_y_init, 50);
title("Histogram of y-Velocities for Electrons in the Semiconductor");
xlabel("y-Velocity (m/s)");
ylabel("Number of Electrons per Bin");
