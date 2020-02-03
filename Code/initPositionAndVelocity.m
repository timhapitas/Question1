function [] = initPositionAndVelocity(randnOrMB, UniformOrFragmented)

global numElectrons;
global x_pos_init;
global y_pos_init;
global V_x_init;
global V_y_init;
global magThermalSpeed;
global boxWidthScaleFactor;
global boxLengthScaleFactor;

if (UniformOrFragmented == "Uniform")
    
    x_pos_init = rand(1, numElectrons)*100*boxWidthScaleFactor;
    y_pos_init = rand(1, numElectrons)*100*boxLengthScaleFactor;

elseif (UniformOrFragmented == "Fragmented")
    
    x_init = rand(1, numElectrons)*100*boxWidthScaleFactor;
    y_init = rand(1, numElectrons)*100*boxLengthScaleFactor;

    while y == x
       
        
    end

end

if(randnOrMB == "randn")

    directionAngles = randn(1, numElectrons)*2*pi;

    V_x_init = magThermalSpeed*cos(directionAngles);
    V_y_init = magThermalSpeed*sin(directionAngles);

elseif(randnOrMB == "MB")   
    
    [V_x_init, V_y_init] = thermalize(numElectrons);
    
else
    
    directionAngles = rand(1, numElectrons)*2*pi;

    V_x_init = magThermalSpeed*cos(directionAngles);
    V_y_init = magThermalSpeed*sin(directionAngles);
    
end

end