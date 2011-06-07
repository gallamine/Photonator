% Rev 4 - switch to split mover/receiver, fix direction vector problem,
% lookup table/interpolation

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

% index from origin, 0 theta (angle between x and z), 0 phi (angle between x and y) along x-axis

clear all
clc

% Program run constants
num_photons = 5e5;
albedo = 0.95;          % Albedo of Maalox (ranges from 0.8 to 0.95)
c = 2.35;
receiver_z = 3.66;  
receiver_x = 0.20;                     % Y position of the receiver (in meters)
receiver_y = 0;                     % Z position of the receiver(in meters)

rec_pos = [0,0; 0,0; 0,0];   % receivers at the same position (on the x/y plane)
rec_aperture = [0.0079, 0.0508, 0.0508]; % rec_aperture = ones(num_rx,1).*0.0508;  % 0.0508 m = 2 inches
rec_fov = [2.27, 0.0785, 0.1745];     % 0.314159 = 18 deg FOV, 2.27 = 130 deg, 0.0785 = 5 deg



photon = zeros(num_photons,8);
photon(:,7) = ones(num_photons,1);  % set weights to one
photon(:,8) = ones(num_photons,1);  % set all active

totaldist = zeros(num_photons,1);   % record total distance traveled by each photon
rec_dist = zeros(num_photons,1);    % total distance photon traveld at point of reception
rec_loc = zeros(num_photons,2);     % location of the received photon on the z,y rec. plane
total_rec_packets = 0;              % Total number of packets to cross detector
total_rec_power = 0;                % Total power to cross detector


% c = 0.5;                            % attenuation coefficient in m^-1
% a = 0.07;                           % absorption coefficient in m^-1
% % harbor water
% c = 2.190;
% a = 0.366;
% [cdf_scatter,angle] = generate_scatter('calc','harbor');

% % coastal water
% c = 0.22 + 0.179;
% a = 0.179;
% [cdf_scatter,angle] = generate_scatter('calc','coastal');
% Clear water
% c = 0.0374 + 0.114;
% a = 0.114;
% [cdf_scatter,angle] = generate_scatter('calc','clear');

% Calculated values
b = c * albedo;
a = c - b;
scattering_events = ceil((c-a)*receiver_z*5)     %five times the scattering attenuation length
if scattering_events < 10
      scattering_events = 10;
end

num_rx = length(rec_pos);

inv_c = 1/c;
inv_b = 1/(c-a);

% Constants and magic numbers
max_uz = 0.99999;
n_water = 1.33;                     % index of refraction of water

[cdf_scatter,angle] = generate_scatter('measured','maalox_alan');

i = 0;

h = waitbar(0.0,'Please wait...','CreateCancelBtn','stop=true; delete(h); clear h');  
set(h,'Name','optional window name');  

scattering_length = 1;
       
% X position of the receiver (in meters)
% receiver_z = 3*inv_c;                     % X position of the receiver (in meters)

% c = scattering_length./receiver_z;
% a = c.*(1-albedo);     

% Photons at origin
photon(:,1) = zeros(num_photons,1);        % 0
photon(:,2) = zeros(num_photons,1);        % 0
photon(:,3) = zeros(num_photons,1);        % 0
% point down z-axis
photon(:,4) = zeros(num_photons,1);        % 0
photon(:,5) = zeros(num_photons,1);        % 0
photon(:,6) = ones(num_photons,1);         % 1

r_avg = 0;

tic;

for j = 1:scattering_events
    rand_array = rand(num_photons,3);    % generate a matrix for each photon with rand propogation, roll, and pitch
    rand_array(:,3) = rand_array(:,3).*2.*pi;   %% Uniformly distributed over 0 -2Pi                                     
    waitbar(j/scattering_events,h,['Scattering event ' num2str(j)]);  
    
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
                   
                % FIX THESE EQUATIONS!!
           
                % z distance to receiver plane
                z_dist_rec_intersection = receiver_z - photon(i,3);
                
                if photon(i,6) ~= 0
                    % y distance to receiver plane
                    y_dist_rec_intersection = z_dist_rec_intersection*photon(i,5)/photon(i,6);      % z * uy/uz = z*tan(theta)*sin(phi)
                    % x distance to receiver plane
                    x_dist_rec_intersection = z_dist_rec_intersection*photon(i,4)/photon(i,6);      % z*ux/uz = z*tan(theta)*cos(phi)
                else
                    disp('how can the photon cross the receiver plane when it"s pointing straight up???');
                end
              

                % euclidian distance to the reciever plane
                %dist_to_rec = sqrt((x_dist_rec_intersection)^2 + (y_dist_rec_intersection)^2 + (z_dist_rec_intersection)^2);
                dist_to_rec = z_dist_rec_intersection / photon(i,6);    % z / mu_z

                
                rec_loc(i,1) = photon(i,1) + x_dist_rec_intersection;   % x-axis location of reception
                rec_loc(i,2) = photon(i,2) + y_dist_rec_intersection;    % y-axis location of reception
                rec_loc(i,3) = acos(photon(i,6));                             % incident angle
                rec_loc(i,4) = photon(i,4);                             % for statistics, should be uniform
                
                total_rec_packets = total_rec_packets + 1;
%                     total_rec_power = total_rec_power + photon(i,7)*exp(-dist_to_rec*a);    
                total_rec_power = total_rec_power + photon(i,7);
                rec_dist(i) = totaldist(i)+ dist_to_rec;
                photon(i,8) = 0;
                
                totaldist(i) = totaldist(i) + dist_to_rec;
                
                

            else        % if the photon didn't move into the detector, reduce its power & move it
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
                
                % update the total distance the photon has traveled 
                totaldist(i) = totaldist(i) + r;
            end
        end
    end
    %r_avg = r_avg/num_photons;
%     figure(1)
%     subplot(3,3,j)
%     scatter(photon(:,1),photon(:,2),50,log10(photon(:,7)),'.')
%     xlim([-15 50]);
%     ylim([-25 25]);
% %     xlabel('x-axis (m)')
% %     xlabel('y-axis (m)')
end





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


[power,ph_cnt,angle_mean,angle_var,power_mean,power_var] = mc_rec_r3(a,rec_loc_final,total_rec_dist,rec_pos,rec_aperture,rec_fov);

power/num_photons               % Total received power / total transmitted photon packets
sd_power = sqrt(power_var).*ph_cnt/num_photons;     % sd_power = sd_photon*rx_photons/total_photons.
    
        
total_time = toc;
minutes_to_run = total_time/60

% figure(8)                   % plot 3D histogram of photons on RX plane
% hist3([rec_loc_final(:,1),rec_loc_final(:,2)],[100 100])
% set(gcf,'renderer','opengl');
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');

% j = 1;
% scattered_times = zeros(total_rec_packets,1);
% for i = 1:num_photons
%     
%    if (photon(i,7) == 0)
%        scattered_times(j) = photon(i,8);
%        j = j + 1;       
%    end
% end

total_rec_power
total_rec_packets
% hist(scattered_times,[0:10]);

close(h)

% figure(10)                              % Plot histogram of time-of-arrival vs. power
% [N,X] = hist(total_rec_dist,100);
% pwr_rx_hist = N.*exp(-a.*X);
% stem(X/(3e8/n_water),pwr_rx_hist/total_rec_power)   
% xlabel('Time of arrival (sec)')
% ylabel('Percentage of power')

% distance_delta = max(X) - min(X);
% time_delta = distance_delta/(3e8/n_water);
% T = mean(X(2:end) - X(1:end-1))/(3e8/n_water); % Effective "sampling rate" of the histogram
% bandwidth = 1/T;                                % Normalized frequency in Hz

%  figure(7)
% freqz(N/total_rec_packets,[1],512,bandwidth)      %% plot a frequency response from the histogram data

% figure(9)                           %% Scatter plot of photons on RX plane
% scatter(rec_loc(:,1),rec_loc(:,2))

% figure(1)
% scatter(photon(:,1),photon(:,2),50,log10(photon(:,7)),'.')
% xlim([-20 150]);
% ylim([-120 120]);
% xlabel('x-axis (m)')
% xlabel('y-axis (m)')
% 
% figure(4)
% scatter3(photon(:,1),photon(:,2),photon(:,3),50,log10(photon(:,7)),'.')

% % Draw the box representing the receiver
% line([receiver_z receiver_z],[receiver_y_min receiver_y_min],[receiver_z_max receiver_z_min],'LineWidth',4)   % |
% line([receiver_z receiver_z],[receiver_y_min receiver_y_max],[receiver_z_max receiver_z_max],'LineWidth',4)   % -
% line([receiver_z receiver_z],[receiver_y_max receiver_y_max],[receiver_z_max receiver_z_min],'LineWidth',4)   %  |
% line([receiver_z receiver_z],[receiver_y_min receiver_y_max],[receiver_z_min receiver_z_min],'LineWidth',4)   % _
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

sprintf('Simulation on DATE with %d photons, %d scattering events.', num_photons, scattering_events)
sprintf('C = %d (1/m), A = %d (1/m).',c,a)
sprintf('Receiver at %d, %d, %d (meters)',receiver_x, receiver_y, receiver_z)
% sprintf('Travel distance delta %d (m). Time of arrival delta %d (sec)', distance_delta, time_delta)
% sprintf('Time delta between histogram bins: %d (sec), %d (Hz)',T,bandwidth)