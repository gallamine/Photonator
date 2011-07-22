function [x,y,ux,uy,uz]  = beamProfile(n,sd,diverg,~)

%sd = 1/e beam radius, E(b) = (1/e)E(0)
muD = cos(diverg);

randVals = rand(n,1);

r = sd*sqrt(-log(1-randVals));
ang = (2*pi).* rand(n,1);

cosPhi = cos(ang);
sinPhi = sin(ang);

x = r.*cosPhi;
y = r.*sinPhi;

uz = 1 - rand(n,1).*(1 - muD); % randomly sample values from [1, muD], ISOTROPIC over [0,acos(muD)] <- IS THIS RIGHT?

sin_uz = sqrt(1-uz.^2);
ux = sin_uz.*cosPhi;
uy = sin_uz.*sinPhi;

% hist3([x y],[100 100])