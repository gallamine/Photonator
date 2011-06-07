% set up an array for the photons, 
% x, y, z, roll, pitch, weight, received
% 0, 0, 0, 0   , 0    , 1     , 1
% 1, 2, 3, 4   , 5    , 6     , 7

% index from origin, 0 theta (angle between x and y), 0 phi (angle between
% x and z) along x-axis

function [total_time,total_rec_power,total_rec_packets,rec_loc_final,total_rec_dist] = ...
    mc_func_r4(num_photons,scattering_events,c,a,receiver_x,...
    cdf_scatter,angle,init_angle,init_angle2,diverg)



photon = zeros(num_photons,7);
photon(:,6) = ones(num_photons,1);  % set weights to one
photon(:,7) = ones(num_photons,1);  % set all active
totaldist = zeros(num_photons,1);   % record total distance traveled by each photon
rec_dist = zeros(num_photons,1);    % total distance photon traveld at point of reception
rec_loc = zeros(num_photons,3);     % location of the received photon on the z,y rec. plane
total_rec_packets = 0;              % Total number of packets to cross detector
total_rec_power = 0;                % Total power to cross detector

prob_of_survival = ((c-a)/c);       % b/c
inv_c = 1/c;
inv_b = 1/(c-a);

% receiver_y_min = receiver_y - aperture/2;
% receiver_y_max = receiver_y + aperture/2;
% receiver_z_min = receiver_z - aperture/2;
% receiver_z_max = receiver_z + aperture/2;



tic;

theta = diverg.*rand(1,num_photons);    % uniform rand. dist over divergence ang.
phi = (2*pi).*rand(1,num_photons);      % uniform rand. dist over azmuthal angles

beta = init_angle;   % direction the transmitter is pointing (zenith)
alpha = init_angle2;     % direction the transmitter is pointing (azmuth)

x = cos(beta).*cos(theta)-sin(beta).*cos(phi).*sin(theta);
y = sin(alpha).*sin(beta).*cos(theta)+(cos(alpha).*sin(phi)+sin(alpha).*cos(beta).*cos(phi)).*sin(theta);
z = cos(alpha).*sin(beta).*cos(theta)+(-sin(alpha).*sin(phi)+cos(alpha).*cos(beta).*cos(phi)).*sin(theta);

photon(:,4) = atan2(sqrt(y.^2+z.^2),x);
photon(:,5) = atan2(y,z);

clear theta phi x y z beta alpha

for j = 1:scattering_events
    rand_array = rand(num_photons,3);    % generate a matrix for each photon with rand propogation, roll, and pitch
    rand_array(:,3) = rand_array(:,3).*2.*pi;   %% Uniformly distributed over 0 -2Pi                                                                             

    
    % iterate over every single photon to calculate new position and
    % whether it was received or not.
    for i = 1:num_photons
        
        % if photon hasn't been received
        if (photon(i,7) == 1)
        
             r = -inv_b*log(rand_array(i,1));     % randomly generate optical path length      

           %Generate scattering angle from beam spread function
%             k = 1;
%             while rand_array(i,2) > cdf_scatter(k)
%                 k = k+1;
%             end
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

            theta = angle(k-1) + (rand_array(i,2) - cdf_scatter(k-1))*(angle(k) - angle(k-1))/(cdf_scatter(k) - cdf_scatter(k-1));   



            %theta = angle(k);

            phi = rand_array(i,3);  % Phi is uniformly distributed over 0-2pi


            % find new position based on PREVIOUS direction of motion. This
            % takes care of the initial condition problem where photons
            % were scattered immediately on the first iteration.
            new_x = r*cos(photon(i,4));
            new_y = r*sin(photon(i,4))*sin(photon(i,5));
            new_z = r*sin(photon(i,4))*cos(photon(i,5));
                
            
%             % if the photon is in the detectors FOV, has the correct
%             % acceptance angle, and is to the left of it
%             if (abs(atan((photon(i,2)-receiver_y)/(cp - photon(i,1)))) <= half_fov ...
%                     && abs(atan((photon(i,3)-receiver_z)/(cp - photon(i,1)))) <= half_fov ...
%                     && abs(photon(i,4)) <= half_fov ...
%                     && photon(i,1) <= receiver_x)

            % if the photon moves past the receiver plane
            if ((photon(i,1) + new_x) >= receiver_x)

           
                % x distance to receiver plane
                x_dist_rec_intersection = receiver_x - photon(i,1);
                % y distance to receiver plane
                y_dist_rec_intersection = x_dist_rec_intersection*tan(photon(i,4))*sin(photon(i,5));
                % z distance to receiver plane
                z_dist_rec_intersection = x_dist_rec_intersection*tan(photon(i,4))*cos(photon(i,5));
              

                % euclidian distance to the reciever plane
                dist_to_rec = sqrt((x_dist_rec_intersection)^2 + (y_dist_rec_intersection)^2 + (z_dist_rec_intersection)^2);

%                 % if the movement vector intersects the receiver plane
%                 %if ((photon(i,1) > receiver_x) && (photon(i,2) - y_dist_rec_intersection) <= receiver_y_max) && (photon(i,2) - y_dist_rec_intersection >= receiver_y_min) && (photon(i,3) - z_dist_rec_intersection <= receiver_z_max) && (photon(i,3) - z_dist_rec_intersection >= receiver_z_min)
%                 if ((dist_to_rec < r)...
%                         && (sqrt(((photon(i,3) + z_dist_rec_intersection) - receiver_z)^2 + ...
%                         ((photon(i,2) + y_dist_rec_intersection) - receiver_y)^2) < rec_radius))
%                     
                    rec_loc(i,1) = photon(i,2) + y_dist_rec_intersection;   % y-axis location of reception
                    rec_loc(i,2) = photon(i,3) + z_dist_rec_intersection;    % z-axis location of reception
                    rec_loc(i,3) = photon(i,4);                             % incident angle
                    
                    total_rec_packets = total_rec_packets + 1;
%                     total_rec_power = total_rec_power + photon(i,6)*exp(-dist_to_rec*a);    
                    total_rec_power = total_rec_power + photon(i,6);
                    rec_dist(i) = totaldist(i)+ dist_to_rec;
                    photon(i,7) = 0;
                    
                    % update the total distance the photon has traveled 
                    totaldist(i) = totaldist(i) + dist_to_rec;
                    
                else % if the photon didn't move into the detector, reduce its power & move it
%                     photon(i,6) = photon(i,6)*exp(-r*a); 
                    %photon(i,6) = photon(i,6)*prob_of_survival;
                    % move to new x position
                    photon(i,1) = photon(i,1) + new_x;                                   
                    % move to new y position
                    photon(i,2) = photon(i,2) + new_y;                                 
                    % move to new z position
                    photon(i,3) = photon(i,3) + new_z;
                    
                    % set new theta angle
                    photon(i,4) = photon(i,4) + theta;
                    if (photon(i,4) > pi)
                        photon(i,4) = photon(i,4) - 2*pi;
                    end
                    % set new phi angle
                    photon(i,5) = photon(i,5) + phi;   
                    
                    % update the total distance the photon has traveled 
                    totaldist(i) = totaldist(i) + r;
                    
%                     if rand > prob_of_survival
%                         photon(i,7) = -1;           % Photon absorbed
%                     end    
                end
% 
%             else                % if the photon isn't in the receiver FOV cone
%                 % set new x position
%                 photon(i,1) = photon(i,1) + new_x;                                   
%                 % set new y position
%                 photon(i,2) = photon(i,2) + new_y;                                 
%                 % set new z position
%                 photon(i,3) = photon(i,3) + new_z;
%                 % update weight
% %                 photon(i,6) = photon(i,6)*exp(-r*a);
%                 %photon(i,6) = photon(i,6)*prob_of_survival;
%                 
%                 % set new theta angle
%                 photon(i,4) = photon(i,4) + theta;
%                 % set new phi angle
%                 photon(i,5) = photon(i,5) + phi;   
%                 
%                 % update the total distance the photon has traveled 
%                 totaldist(i) = totaldist(i) + r;
%                 
% %                 if rand > prob_of_survival
% %                     photon(i,7) = -1;           % Photon absorbed
% %                 end                 
%             end
        end
    end
end

total_time = toc;




rec_loc_final = ones(total_rec_packets,3);
j = 1;
for i = 1:num_photons
    if (photon(i,7) == 0)
       rec_loc_final(j,:) = rec_loc(i,:);
       
       j = j +1;
    end
end


j = 1;
total_rec_dist = zeros(total_rec_packets,1);
total_rec_power = 0;
for i = 1:num_photons
    
   if (photon(i,7) == 0)
       total_rec_dist(j) = rec_dist(i);
       total_rec_power = total_rec_power + exp(-a*rec_dist(i));
       j = j + 1;       
   end
end





