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

dx = zeros(1, numElectrons);
dy = zeros(1, numElectrons);
vx = zeros(1, numElectrons);
vy = zeros(1, numElectrons);
x = zeros(1, numElectrons);
y = zeros(1, numElectrons);

boxWidthScaleFactor = 2e-9; %nm%
boxLengthScaleFactor = 1e-9; %nm%

numElectrons = 1000;
magThermalSpeed = sqrt(pi*C.k_b*C.T/C.m);
maxTimeStep = 1000;

initPositionAndVelocity("rand");

dt = (1/500)*100*boxWidthScaleFactor/magThermalSpeed;

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

v_avg = 0;

for t = 1:maxTimeStep
    
    x_previous = x;
    y_previous = y;
    
    vx_previous = vx;
    vy_previous = vy;
    
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
    
    dx = vx * dt;
    dy = vy * dt;
    
    x = x_previous + dx;
    y = y_previous + dy;
    
    x_oldvals(t,:) = x_previous(electronsToTrack);
    y_oldvals(t,:) = y_previous(electronsToTrack);
    
    v_avg = sqrt(sum((vx.^2) + (vy.^2))/numElectrons);
    T(t) = ((v_avg^2)*C.m)/(pi*C.k_b);
    
    % plot code goes here*
    
    plot(x, y, 'b.');
    axis([0 100*boxWidthScaleFactor 0 100*boxLengthScaleFactor]);
    pause(0.0001);
    
end


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


figure;
plot(t_vec, T, 'r.');
title("Semiconductor Temperature Vs Time");
xlabel("Time (s)");
ylabel("Temperature (K)")
grid on;
