
%%
tic
diverg = 0.001;

thetaD = diverg.*rand(1,2e6);    % uniform rand. dist over divergence ang.
phiD = (2*pi).*rand(1,2e6);      % uniform rand. dist over azmuthal angles

beta = 90.*pi/180;   % direction the transmitter is pointing (zenith)
alpha = 90.*pi/180;     % direction the transmitter is pointing (azmuth)

x = cos(beta).*cos(theta)-sin(beta).*cos(phi).*sin(theta);
y = sin(alpha).*sin(beta).*cos(theta)+(cos(alpha).*sin(phi)+sin(alpha).*cos(beta).*cos(phi)).*sin(theta);
z = cos(alpha).*sin(beta).*cos(theta)+(-sin(alpha).*sin(phi)+cos(alpha).*cos(beta).*cos(phi)).*sin(theta);

thetaF = atan2(sqrt(y.^2+z.^2),x);
phiF = atan2(y,z);

init_angle = 0;     %theta
init_angle2 = 0;     % phi

toc;
% 
% scatter3(x,y,z);
% hold on;
% scatter3([1],[0],[0])

%%
% 
% syms alpha beta phi theta
% 
% Rx_alpha = [1,0,0;0,cos(alpha),sin(alpha);0,-sin(alpha),cos(alpha)];
% Ry_beta = [cos(beta),0,-sin(beta);0,1,0;sin(beta),0,cos(beta)];
% Rx_phi = [1,0,0;0,cos(phi),sin(phi);0,-sin(phi),cos(phi)];
% Ry_theta = [cos(theta),0,-sin(theta);0,1,0;sin(theta),0,cos(theta)];
% init = [1;0;0];
% 
% final = Rx_alpha*Ry_beta*Rx_phi*Ry_theta*init

%%
