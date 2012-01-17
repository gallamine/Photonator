function [x,y,ux,uy,uz]  = beamProfile(n,beamWaist,diverg,~)

% Original beam function was WRONG WRONG WRONG
% Updated 10/31/11 - Start with a Gaussian radius ray bundle distribution
% on the x/y plane. No divergence. Pass this through a THIN LENS ray matrix
% transformation (x1 = x2, theta2 = (-1/f)x1 + theta1) to diverge the beam.
% Choose 1/f to give the desired DIVERGENCE when the ray at the BEAM WAIST
% enters the diverging lens.

randVals = rand(n,1);

invF = -diverg/beamWaist;

radius = beamWaist*sqrt(-log(1-randVals));     % Generate Gaussian radius distribution
divAng = -invF.*radius;                 % Polar divergence angle based on thin lens ray xform matrix
phiAng = (2*pi).* rand(n,1);       	% Uniform distribution azmuth angle   

cosPhi = cos(phiAng);
sinPhi = sin(phiAng);

uz = cos(divAng);
sin_uz = sqrt(1-uz.^2);
ux = sin_uz.*cosPhi;
uy = sin_uz.*sinPhi;

x = radius.*cosPhi;
y = radius.*sinPhi;

% %sd = 1/e beam radius, E(b) = (1/e)E(0)
% muD = cos(diverg);
% 
% randVals = rand(n,1);
% 
% 
% ang = (2*pi).* rand(n,1);
% 
% cosPhi = cos(ang);
% sinPhi = sin(ang);
% 
% x = r.*cosPhi;
% y = r.*sinPhi;
% 
% uz = 1 - rand(n,1).*(1 - muD); % randomly sample values from [1, muD], ISOTROPIC over [0,acos(muD)] <- IS THIS RIGHT?



% hist3([x y],[100 100])