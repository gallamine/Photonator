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
% photon(:,8) == RECEIVED? 1 = No, 0 = Yes

% index from origin, 0 theta (angle between x and z), 0 phi (angle between
% x and y) along x-axis

% V5 - change scattering equations. Change direction of propogation to
% z-axis.

function [total_time,total_rec_power,total_rec_packets,rec_loc_final,total_rec_dist] = ...
    mc_func_r5(num_photons,scattering_events,c,a,receiver_z,...
    cdf_scatter,angle,init_angle,init_angle2,diverg)



photon = zeros(num_photons,8);
photon(:,7) = ones(num_photons,1);  % set weights to one
photon(:,8) = ones(num_photons,1);  % set all active
totaldist = zeros(num_photons,1);   % record total distance traveled by each photon
rec_dist = zeros(num_photons,1);    % total distance photon traveld at point of reception
rec_loc = zeros(num_photons,4);     % location of the received photon on the z,y rec. plane
total_rec_packets = 0;              % Total number of packets to cross detector
total_rec_power = 0;                % Total power to cross detector

prob_of_survival = ((c-a)/c);       % b/c
inv_c = 1/c;
inv_b = 1/(c-a);

attenLen = c*receiver_z;
if (attenLen > 13)                      % exp(-13) ~ 1e-6
    min_power = exp(-1.25*attenLen);    % set the minimum power to be a little less than the power of a balistic photon at the receiver
else
    min_power = 1e-6;                   % minimum power value for photons before they are terminated by rouletting
end


max_uz = 0.99999;

% receiver_y_min = receiver_y - aperture/2;
% receiver_y_max = receiver_y + aperture/2;
% receiver_z_min = receiver_z - aperture/2;
% receiver_z_max = receiver_z + aperture/2;



tic;
% 
% theta = diverg.*rand(1,num_photons);    % uniform rand. dist over divergence ang.
% phi = (2*pi).*rand(1,num_photons);      % uniform rand. dist over azmuthal angles
% 
% beta = init_angle;   % direction the transmitter is pointing (zenith)
% alpha = init_angle2;     % direction the transmitter is pointing (azmuth)

% NEED TO RECALCULATE THIS!!
% x = cos(beta).*cos(theta)-sin(beta).*cos(phi).*sin(theta);
% y = sin(alpha).*sin(beta).*cos(theta)+(cos(alpha).*sin(phi)+sin(alpha).*cos(beta).*cos(phi)).*sin(theta);
% z = cos(alpha).*sin(beta).*cos(theta)+(-sin(alpha).*sin(phi)+cos(alpha).*cos(beta).*cos(phi)).*sin(theta);

% photon(:,4) = atan2(sqrt(y.^2+z.^2),x);
% photon(:,5) = atan2(y,z);

% point down z-axis
photon(:,4) = zeros(num_photons,1);        % 0
photon(:,5) = zeros(num_photons,1);        % 0
photon(:,6) = ones(num_photons,1);         % 1

clear theta phi x y z beta alpha

for j = 1:scattering_events
    rand_array = rand(num_photons,3);    % generate a matrix for each photon with rand propogation, roll, and pitch
    rand_array(:,3) = rand_array(:,3).*2.*pi;   %% Uniformly distributed over (0,2Pi)                                                                             

    
    % iterate over every single photon to calculate new position and
    % whether it was received or not.
    for i = 1:num_photons
        
        % if photon hasn't been received
        if (photon(i,8) == 1)
        
             r = -inv_b*log(rand_array(i,1));     % randomly generate optical path length      

           %Generate scattering angle from beam spread function
            minIndex = 2;
            maxIndex = length(cdf_scatter);
            midIndex = minIndex + ceil((maxIndex - minIndex)/2);

            while rand_array(i,2) ~= cdf_scatter(midIndex) && maxIndex > minIndex

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

%                 current_theta = acos(photon(i,6));      % cos^-1(uz) 
%                 if photon(i,6) ~= 1
%                     if photon(i,5) >= 0
%                         current_phi = acos(photon(i,4)/sqrt(1-photon(i,6)^2));  % cos^-1(ux/sqrt(1-uz^2))
%                     else
%                         current_phi = 2*pi - acos(photon(i,4)/sqrt(1-photon(i,6)^2));  % cos^-1(ux/sqrt(1-uz^2))
%                     end
%                 else
%                     current_phi = rand()*2.*pi;
%                 end
%                 if ~isreal(current_phi)
%                     format long
%                     current_phi 
%                     photon(i,:)
%                     disp('Invalid phi');
%                     current_phi = real(current_phi);
%                 end
           
                if photon(i,6) ~= 0
                    % z distance to receiver plane
                    z_dist_rec_intersection = receiver_z - photon(i,3);
                    % y distance to receiver plane
                    y_dist_rec_intersection = z_dist_rec_intersection*photon(i,5)/photon(i,6);
                    % x distance to receiver plane
                    x_dist_rec_intersection = z_dist_rec_intersection*photon(i,4)/photon(i,6);
                else
                    disp('how can the photon cross the receiver plane when it"s pointing straight up???');
                end

                % euclidian distance to the reciever plane
%                 dist_to_rec = sqrt((x_dist_rec_intersection)^2 + (y_dist_rec_intersection)^2 + (z_dist_rec_intersection)^2);
                dist_to_rec = z_dist_rec_intersection / photon(i,6);    % z / mu_z

                rec_loc(i,1) = photon(i,1) + x_dist_rec_intersection;   % x-axis location of reception
                rec_loc(i,2) = photon(i,2) + y_dist_rec_intersection;    % y-axis location of reception
                rec_loc(i,3) = acos(photon(i,6));                             % incident angle
                rec_loc(i,4) = 0;                             % for statistics, should be uniform

                total_rec_packets = total_rec_packets + 1;
%                     total_rec_power = total_rec_power + photon(i,7)*exp(-dist_to_rec*a);    
                total_rec_power = total_rec_power + photon(i,7);
                rec_dist(i) = totaldist(i)+ dist_to_rec;
                photon(i,8) = 0;
                
%                 if ~isreal(dist_to_rec)
%                     rec_dist(i) 
%                     dist_to_rec
%                     rec_loc(i)
%                     z_dist_rec_intersection
%                     y_dist_rec_intersection
%                     x_dist_rec_intersection
%                     current_phi
%                     current_theta
%                     error('Invalid rec dist');
%                 end

                % update the total distance the photon has traveled 
                totaldist(i) = totaldist(i) + dist_to_rec;

            else % if the photon didn't move into the detector, reduce its power & move it
                %photon(i,7) = photon(i,7)*exp(-r*a); 
                %photon(i,7) = photon(i,7)*prob_of_survival;
                % move to new x position
                photon(i,1) = photon(i,1) + x_step;                                   
                % move to new y position
                photon(i,2) = photon(i,2) + y_step;                                 
                % move to new z position
                photon(i,3) = photon(i,3) + z_step;


                if abs(photon(i,6)) > max_uz         % if uz ~ 1
                    %photon(i,:)
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

%                 if ((photon(i,6) > 1) || (photon(i,5) > 1) || (photon(i,4) > 1))
%                     j
%                      theta
%                      phi
%                      old_ux
%                      old_uy
%                      old_uz
%                      photon(i,:)
%                     error('Invalid cosine');
%                 end
                
%                 if (photon(i,6)^2 + photon(i,5)^2 + photon(i,4)^2) < 0.90
%                     j
%                     diff = 1- (photon(i,6)^2 + photon(i,5)^2 + photon(i,4)^2)
%                     theta
%                     phi
%                     ph = photon(i,:)
%                     error('Invalid cosine vecotrs');
%                 end

                % update the total distance the photon has traveled 
                totaldist(i) = totaldist(i) + r;

%                     if rand > prob_of_survival
%                         photon(i,8) = -1;           % Photon absorbed
%                     end    
            end
        end
    end
end

total_time = toc;



% Do this after the loop above, so we can allocate the array once (instead
% of growing dynamically as we receive individual photons - SLOW)
rec_loc_final = ones(total_rec_packets,4);
j = 1;
for i = 1:num_photons
    if (photon(i,8) == 0)
       rec_loc_final(j,:) = rec_loc(i,:);
       
       j = j +1;
    end
end


j = 1;
total_rec_dist = zeros(total_rec_packets,1);
total_rec_power = 0;
for i = 1:num_photons
    
   if (photon(i,8) == 0)
       total_rec_dist(j) = rec_dist(i);
       total_rec_power = total_rec_power + exp(-a*rec_dist(i));
       j = j + 1;       
   end
end





