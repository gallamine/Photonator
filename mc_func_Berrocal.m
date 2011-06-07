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
    mc_func_Berrocal(num_photons,scattering_events,c,a,receiver_z,...
    cdf_scatter,angle,init_angle,init_angle2,diverg)

useLimits = 'true';

if strcmp(useLimits,'true')
    xLimMax = 0.005;
    xLimMin = -0.005;
    yLimMax = 0.005;
    yLimMin = -0.005;
    zLimMax = 0.01;
    zLimMin = -0.005;
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
inv_b = 1/(c-a);

min_power = 1e-4;                   % minimum power value for photons before they are terminated by rouletting

max_uz = 0.99999;

tic;
% 
% theta = diverg.*rand(1,num_photons);    % uniform rand. dist over divergence ang.
% phi = (2*pi).*rand(1,num_photons);      % uniform rand. dist over azmuthal angles
% 
% beta = init_angle;   % direction the transmitter is pointing (zenith)
% alpha = init_angle2;     % direction the transmitter is pointing (azmuth)

% point down z-axis
[x,y] = startingDistribution(num_photons);         % Assuming a square end size
photon(:,1) = xLimMin;                              % x position
photon(:,2) = y.*yLimMax;                           % y position
photon(:,3) = x.*xLimMax;                           % z position

photon(:,4) = ones(num_photons,1);         % z - 1
photonsRemaining = num_photons;             % count down for each photon received/terminated

clear theta phi  z beta alpha

% for j = 1:scattering_events
while photonsRemaining > 0                      % which is faster? create random values on the fly??  <- check
    rand_array = rand(num_photons,3);           % generate a matrix for each photon with rand propogation, roll, and pitch
    rand_array(:,3) = rand_array(:,3).*2.*pi;   % Uniformly distributed over (0,2Pi)                                                                             
    
    % iterate over every single photon to calculate new position and
    % whether it was received or not.
    for i = 1:num_photons
        
        % if photon hasn't been received
        if (photon(i,8) == 1)
        
             r = -inv_c*log(rand_array(i,1));     % randomly generate optical path length      

            %Generate scattering angle from beam spread function
            minIndex = 2;
            maxIndex = length(cdf_scatter);
            midIndex = minIndex + ceil((maxIndex - minIndex)/2);

            while rand_array(i,2) ~= cdf_scatter(midIndex) && maxIndex >= minIndex        % changed > to >=      
    
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
                

            % if the photon moves past the receiver plane
            if ((photon(i,3) + z_step) >= receiver_z)
          
                if photon(i,6) ~= 0                         % If the photon has a z-component, mu_z != 0
                    % z distance to receiver plane
                    z_dist_rec_intersection = receiver_z - photon(i,3);     % Zrec - Zphoton
                    % y distance to receiver plane
                    y_dist_rec_intersection = z_dist_rec_intersection*photon(i,5)/photon(i,6);
                    % x distance to receiver plane
                    x_dist_rec_intersection = z_dist_rec_intersection*photon(i,4)/photon(i,6);
                else
                    disp('how can the photon cross the receiver plane when it"s pointing straight up??? Z distance step:');
                    disp(z_step);
                    disp('Z position:');
                    disp(photon(i,3));
                end

                % euclidian distance to the reciever plane
                dist_to_rec = z_dist_rec_intersection / photon(i,6);    % z / mu_z

                rec_loc(i,1) = photon(i,1) + x_dist_rec_intersection;   % x-axis location of reception
                rec_loc(i,2) = photon(i,2) + y_dist_rec_intersection;    % y-axis location of reception
                rec_loc(i,3) = photon(i,4);                             % incident angle (mu_x)
                rec_loc(i,4) = photon(i,5);                             % for statistics, should be uniform (mu_y)
                rec_loc(i,5) = photon(i,6);                             % for statistics, should be uniform (mu_z)

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
                
                if ((photon(i,1) > xLimMax) || (photon(i,1) < xLimMin) || (photon(i,2) > yLimMax) || (photon(i,2) < yLimMin) || (photon(i,3) < zLimMin))
                    photon(i,8) = -1;                           % mark as terminated
                    photonsRemaining = photonsRemaining - 1;    % decrement outstanding photons
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
                        %disp('mu_z less than 1')
                    end
                end
                % update the total distance the photon has traveled 
                totaldist(i) = totaldist(i) + r;

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




