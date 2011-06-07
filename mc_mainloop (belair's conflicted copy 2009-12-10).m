% set up an array for the photons, 
% x, y, z, roll, pitch, weight
% 0, 0, 0, 0   , 0    , 1
% 1, 2, 3, 4   , 5    , 6

% index from origin, 0 theta (angle between x and y), 0 phi (angle between
% x and z) along x-axis
clear all
clc
[cdf_scatter,angle] = generate_scatter(0.924);

num_photons = 1e6;
scattering_events = 3;

photon = zeros(num_photons,7);
photon(:,6) = ones(num_photons,1);  % set weights to one
photon(:,7) = ones(num_photons,1);  % set all active
totaldist = zeros(num_photons,1);   % record total distance traveled by each photon
rec_dist = zeros(num_photons,1);    % total distance photon traveld at point of reception
total_rec_packets = 0;              % Total number of packets to cross detector
total_rec_power = 0;                % Total power to cross detector


c = 0.8;                            % attenuation coefficient in m^-1
a = 0.09;                           % absorption coefficient in m^-1
inv_c = 1/c;

i = 0;

h = waitbar(0.0,'Please wait...','CreateCancelBtn','stop=true; delete(h); clear h');  
set(h,'Name','optional window name');  

receiver_x = 5;                     % X position of the receiver (in meters)
receiver_y = 0;                     % Y position of the receiver (in meters)
receiver_z = 0;                     % Z position of the receiver(in meters)
aperture = 0.2;                       % Diameter of aperture (in meters)
receiver_y_min = receiver_y - aperture/2;
receiver_y_max = receiver_y + aperture/2;
receiver_z_min = receiver_z - aperture/2;
receiver_z_max = receiver_z + aperture/2;
fov = pi/4;                           % Field of view (in radians)

half_fov = fov/2;

cp = receiver_x + (aperture/2)/tan(half_fov);    % point of receiver FOV 'cone' on the x-axis - different from focal point

init_angle = atan(receiver_y/receiver_x)       % Point transmitter at receiver
photon(:,4) = init_angle;
init_angle2 = atan(receiver_z/receiver_x)       % Point transmitter at receiver
photon(:,5) = init_angle2;


tic;
for j = 1:scattering_events
    rand_array = rand(num_photons,3);    % generate a matrix for each photon with rand propogation, roll, and pitch
                                            
    waitbar(j/scattering_events,h,['Scattering event ' num2str(j)]);  
    
    for i = 1:num_photons
        
        % if photon hasn't been received
        if (photon(i,7) == 1)
        
             r = -inv_c*log(rand_array(i,1));     % randomly generate optical path length      

           %Generate scattering angle from beam spread function
           theta = angle(find(rand_array(i,2) <= cdf_scatter,1));
           phi = angle(find(rand_array(i,3) <= cdf_scatter,1));

    %         % Generate generic  spread
    %         theta = 0.25*pi*rand - 0.125*pi;
    %         phi = 0.25*pi*rand - 0.125*pi;

%             %No scattering
%             theta = 0;
%             phi = 0;


            % calculate new position of photon based on scattering length and
            % propagation direction
            new_x = r*cos(photon(i,5) + phi)*cos(photon(i,4) + theta);
            new_y = r*cos(photon(i,5) + phi)*sin(photon(i,4) + theta);
            new_z = r*cos(photon(i,4) + theta)*sin(photon(i,5) + phi);

            % set new theta angle
            photon(i,4) = photon(i,4) + theta;
            % set new phi angle
            photon(i,5) = photon(i,5) + phi;   
                
       
            % if the photon is in the detectors FOV, has the correct
            % acceptance angle, and is to the left of it
            if (abs(atan((photon(i,2)-receiver_y)/(cp - photon(i,1)))) < half_fov ...
                    && abs(atan((photon(i,3)-receiver_z)/(cp - photon(i,1)))) < half_fov ...
                    && abs(photon(i,4)) <= half_fov ...
                    && abs(photon(i,5)) <= half_fov ...
                    && photon(i,1) <= receiver_x)

% 
%             % if the photon has the correct acceptance angle, and is to
%             % the left of the receiver
%             if (abs(photon(i,4)) <= half_fov ...
%               && abs(photon(i,5)) <= half_fov ...
%               && photon(i,1) <= receiver_x)
%                 
%             if (photon(i,1) <= receiver_x)
            
                % x distance to receiver plane
                x_dist_rec_intersection = receiver_x - photon(i,1);
                % y distance to receiver plane
                y_dist_rec_intersection = x_dist_rec_intersection*tan(photon(i,4));
                % z distance to receiver plane
                z_dist_rec_intersection = x_dist_rec_intersection*tan(photon(i,5));
              

                % euclidian distance to the reciever plane
                dist_to_rec = sqrt((x_dist_rec_intersection)^2 + (y_dist_rec_intersection)^2 + (z_dist_rec_intersection)^2);

                % if the movement vector intersects the receiver plane
                %if ((photon(i,1) > receiver_x) && (photon(i,2) - y_dist_rec_intersection) <= receiver_y_max) && (photon(i,2) - y_dist_rec_intersection >= receiver_y_min) && (photon(i,3) - z_dist_rec_intersection <= receiver_z_max) && (photon(i,3) - z_dist_rec_intersection >= receiver_z_min)
                if ((dist_to_rec < r)...
                        && ((photon(i,2) + y_dist_rec_intersection) <= receiver_y_max)...   % can probably remove these conditions, since we already check the acceptance angle
                        && ((photon(i,2) + y_dist_rec_intersection) >= receiver_y_min)...
                        && ((photon(i,3) + z_dist_rec_intersection) <= receiver_z_max)...
                        && ((photon(i,3) + z_dist_rec_intersection) >= receiver_z_min))
                    % sqrt(((photon(i,2) + y_dist_rec_intersection) - receiver_y)^2 ...
                    % + ((photon(i,3) + z_dist_rec_intersection) -
                    % receiver_z)^2) < aperture/2
                    
                    total_rec_packets = total_rec_packets + 1;
                    total_rec_power = total_rec_power + photon(i,6)*exp(-dist_to_rec*a);                    
                    rec_dist(i) = totaldist(i)+ dist_to_rec;
                    photon(i,7) = 0;
                    
                    % update the total distance the photon has traveled 
                    totaldist(i) = totaldist(i) + dist_to_rec;
                    
                else % if the photon didn't move into the detector, reduce its power & move it
                    photon(i,6) = photon(i,6)*exp(-r*a); 
                    % move to new x position
                    photon(i,1) = photon(i,1) + new_x;                                   
                    % move to new y position
                    photon(i,2) = photon(i,2) + new_y;                                 
                    % move to new z position
                    photon(i,3) = photon(i,3) + new_z;
                    
                    % update the total distance the photon has traveled 
                    totaldist(i) = totaldist(i) + r;
                end

            else                % if the photon isn't in the receiver FOV cone
                % set new x position
                photon(i,1) = photon(i,1) + new_x;                                   
                % set new y position
                photon(i,2) = photon(i,2) + new_y;                                 
                % set new z position
                photon(i,3) = photon(i,3) + new_z;
                % update weight
                photon(i,6) = photon(i,6)*exp(-r*a);
                
                % update the total distance the photon has traveled 
                totaldist(i) = totaldist(i) + r;
            end
        end
    end
end

total_time = toc;
minutes_to_run = total_time/60

j = 1;
total_rec_dist = zeros(total_rec_packets,1);
for i = 1:num_photons
    
   if (photon(i,7) == 0)
       total_rec_dist(j) = totaldist(i);
       j = j + 1;       
   end
end

total_rec_power
total_rec_packets

close(h)

figure(10)
[N,X] = hist(total_rec_dist,100);
plot(X/3e8,N/total_rec_packets)

% figure(1)
% scatter(photon(:,1),photon(:,2),50,log10(photon(:,6)),'.')
% xlim([-20 150]);
% ylim([-120 120]);
% xlabel('x-axis (m)')
% xlabel('y-axis (m)')
% 
% figure(4)
% scatter3(photon(:,1),photon(:,2),photon(:,3),50,log10(photon(:,6)),'.')
% 
% % Draw the box representing the receiver
% line([receiver_x receiver_x],[receiver_y_min receiver_y_min],[receiver_z_max receiver_z_min],'LineWidth',4)   % |
% line([receiver_x receiver_x],[receiver_y_min receiver_y_max],[receiver_z_max receiver_z_max],'LineWidth',4)   % -
% line([receiver_x receiver_x],[receiver_y_max receiver_y_max],[receiver_z_max receiver_z_min],'LineWidth',4)   %  |
% line([receiver_x receiver_x],[receiver_y_min receiver_y_max],[receiver_z_min receiver_z_min],'LineWidth',4)   % _
% xlabel('x-axis (m)')
% ylabel('y-axis (m)')
% zlabel('z-axis (m)')

% figure(3);
% hist(photon(:,1),max(photon(:,1)));


% figure(2)
% hist(totaldist,20);
% figure(3)
% hist(photon(:,4),20)
findfigs