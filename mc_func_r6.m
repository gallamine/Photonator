% set up an array for the photons, 
% x, y, z, mu_x, mu_y, mu_z, weight, received
% 0, 0, 0, 0   , 0   , 1   , 1     , 1  - initial values
% 1, 2, 3, 4   , 5   , 6   , 7     , 8  - Position 
% photon(:,1) == X POSITION
% photon(:,2) == Y POSITION
% photon(:,3) == Z POSITION
% photon(:,4) == ux
% photon(:,5) == uy
% photon(:,6) == uz
% photon(:,7) == WEIGHT
% photon(:,8) == RECEIVED? 1 = No, 0 = Yes (detected), -1 = terminated

% index from origin, 0 theta (angle between x and z), 0 phi (angle between
% x and y) along x-axis

% V5 - change scattering equations. Change direction of propogation to
% z-axis.
% V6 - Lots of changes:
%   - Added size limits to the simulation. Photons terminate when
%   encountering boundaris
%   - Added unlimited scattering. Scattering is dependant on minimum weight
%   - added rouletting to remove photons of low weight
%   - added window effects to the receiver
%   - changed events to be based on 'c', not 'b'
%   - changed how weights are calculated/updated

function [total_time,total_rec_power,total_rec_packets,rec_loc_final,total_rec_dist,rec_weights] = ...
    mc_func_r6(num_photons,scattering_events,c,a,receiver_z,...
    cdf_scatter,angle,init_angle,init_angle2,beamDiverg,beamWidth,wallAbsorption)

useLimits = 'false';
reflection = 0;

rxPlaneLimits = 1;             % Use dimention limits on the receiver plane. This reduces the size of the dataset
rxXLimMax = 3;          % 3 meters
rxXLimMin = -3;
rxYLimMax = 3;
rxYLimMin = -3; 

nAir = 1;
nWater = 1.33;
surfCritAngCos = sqrt(1 - (nAir/nWater)^2);             % sine T1 = nair/nwater sine T3

sizeMult = 1;

zLimMin = 0;

if strcmp(useLimits,'true')
%     xLimMax = 0.61*sizeMult;
%     xLimMin = -.61*sizeMult;
%     yLimMax = 0.61*sizeMult;
%     yLimMin = -.61*sizeMult;
%     zLimMax = receiver_z;
%     zLimMin = 0;

    % Actual tank dimensions from exp in 7_11
    xLimMax = 0.61*sizeMult;
    xLimMin = -.61*sizeMult;
    yLimMax = 0.1905*sizeMult;
    yLimMin = -.5842*sizeMult;
    zLimMax = receiver_z;
    
    
end

photon = zeros(num_photons,8);
photon(:,7) = ones(num_photons,1);  % set weights to one
photon(:,8) = ones(num_photons,1);  % set all active
totaldist = zeros(num_photons,1);   % record total distance traveled by each photon
rec_dist = zeros(num_photons,1);    % total distance photon traveld at point of reception
rec_loc = zeros(num_photons,4);     % location of the received photon on the z,y rec. plane
total_rec_packets = 0;              % Total number of packets to cross detector
total_rec_power = 0;                % Total power to cross detector

prob_of_survival = ((c-a)/c);       % b/c
rouletteConst = 10;                 % 1/rouletteConst determines the probability of a low power photon surviving
rouletteConst_inv = 1/rouletteConst;
inv_c = 1/c;


if prob_of_survival >= 0.90
    min_power = 1e-4;                   % minimum power value for photons before they are terminated by rouletting
    %min_power = 0.5
elseif prob_of_survival >= 0.8299
    min_power = 1e-5;  
elseif prob_of_survival >= 0.70
    min_power = 1e-5;  
else
    min_power = 1e-7;
end
max_uz = 1-1e-12;

tic;
% 
% theta = diverg.*rand(1,num_photons);    % uniform rand. dist over divergence ang.
% phi = (2*pi).*rand(1,num_photons);      % uniform rand. dist over azmuthal angles
% 
% beta = init_angle;   % direction the transmitter is pointing (zenith)
% alpha = init_angle2;     % direction the transmitter is pointing (azmuth)

[photon(:,1),photon(:,2),photon(:,4),photon(:,5),photon(:,6)] = beamProfile(num_photons,beamWidth,beamDiverg,'gaussian');

% % point down z-axis
% photon(:,4) = zeros(num_photons,1);        % x - 0
% photon(:,5) = zeros(num_photons,1);        % y - 0
% photon(:,6) = ones(num_photons,1);         % z - 1
photonsRemaining = num_photons;             % count down for each photon received/terminated

clear theta phi x y z beta alpha

% for j = 1:scattering_events
while photonsRemaining > 0                      % which is faster? create random values on the fly??  <- check
    rand_array = rand(num_photons,3);           % generate a matrix for each photon with rand propogation, roll, and pitch
    rand_array(:,3) = rand_array(:,3).*2.*pi;   % Uniformly distributed over (0,2Pi)                                                                             
    
    % iterate over every single photon to calculate new position and
    % whether it was received or not.
    for i = 1:num_photons
        
        % if the photon is still active (1 - active, 0 - received, -1 - terminated
        if (photon(i,8) == 1)
        
            r = -inv_c*log(rand_array(i,1));     % randomly generate optical path length      

            %Generate scattering angle from beam spread function
            minIndex = 2;
            maxIndex = length(cdf_scatter);
            midIndex = minIndex + ceil((maxIndex - minIndex)/2);

            while rand_array(i,2) ~= cdf_scatter(midIndex) && maxIndex >= minIndex          % Changed > to >=

                midIndex = minIndex + ceil((maxIndex - minIndex)/2);
                if rand_array(i,2) > cdf_scatter(midIndex)
                    minIndex = midIndex + 1;
                elseif rand_array(i,2) < cdf_scatter(midIndex)
                    maxIndex = midIndex - 1;
                end

            end
            midIndex = minIndex + ceil((maxIndex - minIndex)/2);
            k = midIndex;

            theta = angle(k-1) + (rand_array(i,2) - cdf_scatter(k-1))*(angle(k) - angle(k-1))/(cdf_scatter(k) - cdf_scatter(k-1));   % Linear interpolation from "angle" vector
            phi = rand_array(i,3);  % Phi is uniformly distributed over 0-2pi


            % find new position increment based on PREVIOUS direction vectors (ux,uy,uz). This
            % takes care of the initial condition problem where photons
            % were scattered immediately on the first iteration.
            x_step = r*photon(i,4);  % x step size
            y_step = r*photon(i,5);  % y step size
            z_step = r*photon(i,6);  % z step size
                

            % if the photon moves past the receiver plane <- photons that
            % should be reflected are not counted properly (i.e. they
            % aren't reflected before being received)
            if ((photon(i,3) + z_step) >= receiver_z)
          
                if photon(i,6) ~= 0                         % If the photon has a z-component, mu_z != 0
                    % z distance to receiver plane
                    z_dist_rec_intersection = receiver_z - photon(i,3);     % Zrec - Zphoton
                    % y distance to receiver plane
                    y_dist_rec_intersection = z_dist_rec_intersection*photon(i,5)/photon(i,6);
                    % x distance to receiver plane
                    x_dist_rec_intersection = z_dist_rec_intersection*photon(i,4)/photon(i,6);
                else
                    disp('how can the photon cross the receiver plane when it"s pointing straight up???');
                end

                % euclidian distance to the reciever plane
                dist_to_rec = z_dist_rec_intersection / photon(i,6);    % z / mu_z
                
                if rxPlaneLimits == 1                   % If we're limiting the dimensions of the receiver plane (to reduce the dataset)
                    if ((photon(i,1) + x_dist_rec_intersection) > rxXLimMax || (photon(i,1) + x_dist_rec_intersection) < rxXLimMin ...
                            || photon(i,2) + y_dist_rec_intersection > rxYLimMax || photon(i,2) + y_dist_rec_intersection < rxYLimMin)  % Photon exceeds the limits
                        photon(i,8) = -1;                           % mark as terminated
                        photonsRemaining = photonsRemaining - 1;    % decrement outstanding photons
                        continue;                                   % Continue on to the next photon, skipping the code below
                    end
                end

                rec_loc(i,1) = photon(i,1) + x_dist_rec_intersection;   % x-axis location of reception
                rec_loc(i,2) = photon(i,2) + y_dist_rec_intersection;    % y-axis location of reception
                rec_loc(i,3) = photon(i,4);                             % for statistics, should be uniform (mu_x)
                rec_loc(i,4) = photon(i,5);                             % for statistics, should be uniform (mu_y)
                rec_loc(i,5) = photon(i,6);                             % incident angle, mu_z

                total_rec_packets = total_rec_packets + 1;
                total_rec_power = total_rec_power + photon(i,7);        % total power at the receiver plane (sum of received photons)
                
                rec_dist(i) = totaldist(i)+ dist_to_rec;                % individual photon's distance traveled.
                photon(i,8) = 0;                                        % mark photon as received
                photonsRemaining = photonsRemaining - 1;                % decrement number of outstanding photons
                
                % update the total distance the photon has traveled 
                totaldist(i) = totaldist(i) + dist_to_rec;

            else % if the photon didn't move into the detector, reduce its power & move it
                
                photon(i,1) = photon(i,1) + x_step;         % move to new x position                            
                photon(i,2) = photon(i,2) + y_step;         % move to new y position                        
                photon(i,3) = photon(i,3) + z_step;         % move to new z position
                
                 % update the total distance the photon has traveled 
                totaldist(i) = totaldist(i) + r;
                
                if (reflection == 1)                        % Check to make sure we're using reflection in our simulation (perfect reflectors)
                    % While any of the boundaries are exceeded, reflect
                    % around that boundary. Keep doing this till the photon
                    % is contained in the boundaries.
                    while ((photon(i,1) > xLimMax) || (photon(i,1) < xLimMin)  || (photon(i,2) > yLimMax) || (photon(i,2) < yLimMin))
                        if (photon(i,2) > yLimMax)                    % If the photon leaves the top of the volume, reflect downwards
                            
                           if (photon(i,5) < surfCritAngCos)    % TOTAL INTERNAL REFLECTION
                                photon(i,5) = -1*photon(i,5);               % reflect the light beam by flipping the sign of the mu_y vector
                                photon(i,2) = yLimMax - (photon(i,2) - yLimMax);      % New y_pos is Y-Max - (Y_pos - Ymax). Remove the part that "sticks out" of the top of the surface.
                           else            % OTHERWISE, Lose power out the top                    
                               % Calculate fresnel loss
                               cosExitAng = sqrt(1 - (nWater/nAir)^2*(1-photon(i,5)^2));
                               rp = (photon(i,5) - nWater*cosExitAng)/(photon(i,5) + nWater*cosExitAng); 	% Thanks to -> Small Monte Carlo by Scott Prahl (http://omlc.ogi.edu)
                               rs = (cosExitAng - nWater*photon(i,5))/(cosExitAng + nWater*photon(i,5));
                               R = (rp^2 + rs^2)/2;                                                     % unpolarized reflection coefficient
                               photon(i,7) = photon(i,7)*R;                 % Reduce the photon weight        

                                photon(i,5) = -1*photon(i,5);               % reflect the light beam by flipping the sign of the mu_y vector
                                % X and Z position stay the same
                                photon(i,2) = yLimMax - (photon(i,2) - yLimMax);      % New y_pos is Y-Max - (Y_pos - Ymax). Remove the part that "sticks out" of the top of the surface.
                           end
                       elseif (photon(i,2) < yLimMin)                  % Leaves the bottom of the container
                            photon(i,5) = -1*photon(i,5);               % reflect the light beam by flipping the sign of the mu_y vector
                            photon(i,2) = yLimMin - (photon(i,2) - yLimMin);
                            photon(i,7) = photon(i,7)*wallAbsorption;
%                             photon(i,8) = -1;                           % mark as terminated
%                             photonsRemaining = photonsRemaining - 1;    % decrement outstanding photons
%                             continue;
                        end
                        
                        if (photon(i,1) < xLimMin)                   % Leaves the side of the container
                            photon(i,4) = -1*photon(i,4);               % reflect the light beam by flipping the sign of the mu_x vector
                            photon(i,1) = xLimMin - (photon(i,1) - xLimMin);
                            photon(i,7) = photon(i,7)*wallAbsorption;
%                             photon(i,8) = -1;                           % mark as terminated
%                             photonsRemaining = photonsRemaining - 1;    % decrement outstanding photons
%                             continue;
                        elseif (photon(i,1) > xLimMax)                  % Leaves the side of the container
                            photon(i,4) = -1*photon(i,4);               % reflect the light beam by flipping the sign of the mu_x vector
                            photon(i,1) = xLimMax - (photon(i,1) - xLimMax);
                            photon(i,7) = photon(i,7)*wallAbsorption;
%                             photon(i,8) = -1;                           % mark as terminated
%                             photonsRemaining = photonsRemaining - 1;    % decrement outstanding photons
%                             continue;
                        end
                    end
                end                                    
                    
                % If the photon collides with the back wall, terminate
                if (photon(i,3) < zLimMin)
                    photon(i,8) = -1;                           % mark as terminated
                    photonsRemaining = photonsRemaining - 1;    % decrement outstanding photons
                    continue;
                else 	% if the photon is still in the boundaries                      
               
                    photon(i,7) = photon(i,7)*prob_of_survival;     % reduce weight
                    if  photon(i,7) < min_power
                        if rand() > (rouletteConst_inv)                 % if the photon is randomly chosen to be terminated ...
                            photon(i,8) = -1;                               % -1 indicates a photon was terminated, but not received
                            photonsRemaining = photonsRemaining - 1;        % decrement outstanding photons
                        else                                            % ... otherwise the photon gets the energey of terminated photons
                            photon(i,7) = photon(i,7)*rouletteConst;    % shift power of terminated photons to new photon
                        end
                    end
                    
                    if abs(photon(i,6)) > max_uz         % if uz ~ 1
                        photon(i,4) = sin(theta)*cos(phi);
                        photon(i,5) = sin(theta)*sin(phi);
                        photon(i,6) = sign(photon(i,6))*cos(theta);
                        %disp('mu_z near 1')
                    else
                        sqrt_uz = sqrt(1 - photon(i,6)^2);
                        old_ux = photon(i,4);
                        old_uy = photon(i,5);
                        old_uz = photon(i,6);
                        photon(i,4) = (sin(theta)/sqrt_uz)*(old_ux*old_uz*cos(phi) - old_uy*sin(phi)) + old_ux*cos(theta);   % ux
                        photon(i,5) = (sin(theta)/sqrt_uz)*(old_uy*old_uz*cos(phi) + old_ux*sin(phi)) + old_uy*cos(theta);   % uy
                        photon(i,6) = (-sin(theta)*cos(phi))*sqrt_uz + old_uz*cos(theta);                                    % uz
                    end
                    
                    % Normalize the pointing vectors -> ux^2 + uy^2 + uz^2 = 1^2
                    if abs(1 - (photon(i,4)^2 + photon(i,5)^2 + photon(i,6)^2)) > 1e-11
                        normLength = sqrt(photon(i,4)^2 + photon(i,5)^2 + photon(i,6)^2);
                        photon(i,4) = photon(i,4) / normLength;
                        photon(i,5) = photon(i,5) / normLength;
                        photon(i,6) = photon(i,6) / normLength;
                        %disp('Vector normalization wrong!');
                    end                    
                end
            end
        end
    end
end

total_time = toc;

% Do this after the loop above, so we can allocate the array once (instead
% of growing dynamically as we receive individual photons - SLOW)
rec_loc_final = ones(total_rec_packets,5);
j = 1;
for i = 1:num_photons                           % iterate over all photons
    if (photon(i,8) == 0)                       % if the photon was received
       rec_loc_final(j,:) = rec_loc(i,:);       % record the receive location and angles
       j = j +1;                                % increment the number of received photons
    end
end

if ((j-1) ~= total_rec_packets)
    disp('Error! Total received photons doesnt equal j iterator. ');
    disp(sprintf('j = %d and total_rec_packets = %d',j, total_rec_packets));
end

j = 1;
total_rec_dist = zeros(total_rec_packets,1);
rec_weights = zeros(total_rec_packets,1);
total_rec_power = 0;
for i = 1:num_photons    
   if (photon(i,8) == 0)
       total_rec_dist(j) = rec_dist(i);         % store the distance traveled by each photon
       rec_weights(j) = photon(i,7);            % store the weights of received photons
       j = j + 1;       
   end
end

if ((j-1) ~= total_rec_packets)
    disp('Error! Total received distances doesnt equal j iterator');
    disp(sprintf('j = %d and total_rec_packets = %d',j, total_rec_packets));
end




