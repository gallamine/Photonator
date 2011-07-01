function [x,y,uz]  = beamProfile(n,sd,diverg,~)

randVals = rand(n,1);

r = sqrt(-log(randVals).*(sd^2)/2);
ang = (2*pi).* rand(n,1);
divAng = diverg .* rand(n,1);

x = r.*cos(ang);
y = r.*sin(ang);

uz = cos(divAng);       % Return the projecton onto the transmission axis

% hist3([x y],[100 100])