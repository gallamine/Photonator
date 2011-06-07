% VERSION LOG
% Rev 1 - first creation, record power and count for each receiver
% Rev 2 - add statistics, mean and variance for power, angle(s)
% Rev 3 - changed receiver plane to x/y instead of y/z - photons move along
% the z-axis now

function [power,ph_cnt,angle_mean,angle_var,power_mean,power_var] = mc_rec_r3(a,rec_loc_final,total_rec_dist,rec_pos,rec_aperture,rec_fov)

num_rx = length(rec_pos);   % total number of receivers
num_photons = length(rec_loc_final);    % total number of received photons

ph_cnt = zeros(1,num_rx);
power = zeros(1,num_rx);

angle_mean = zeros(1,num_rx);   % mean of the incident angle for each Rx
angle_var = zeros(1,num_rx);    % variance of incident angle for each Rx
power_mean = zeros(1,num_rx);    % Mean distance traveled for each Rx
power_var = zeros(1,num_rx);     % Variance of dist traveled for each Rx

dWindow = 0.00635;              % 0.25 inches
dAir = 0.127;                   % 5 inches
nWater = 1.33;
nWindow = 1.585;
nAir = 1;

nWaterWindow = nWater/nWindow;
nWindowAir = nWindow/nAir;
   
for i = 1:num_photons
   ph_x = rec_loc_final(i,1);
   ph_y = rec_loc_final(i,2);
   mu_x = rec_loc_final(i,3);
   mu_y = rec_loc_final(i,4); 
   mu_z = rec_loc_final(i,5); 
   
   
   for j = 1:num_rx
       rx_x = rec_pos(j,1);
       rx_y = rec_pos(j,2);
       radius = rec_aperture(j)/2;
       angle_sqravg = 0;
       power_sqravg = 0;   
       
       % Find new position from window
       sinT = sqrt(1 - mu_z^2);
       sinT1 =  nWaterWindow*sinT;          % 1.33/1.585 * sqrt(1 - uz^2)
       mu_z1 = sqrt(1-sinT1^2);
       mu_x1 = nWaterWindow * mu_x;
       mu_y1 = nWaterWindow * mu_y;
       delY1 = dWindow * mu_y1/mu_z1;
       delX1 = dWindow * mu_x1/mu_z1;
       
       sinT2 = nWater*sinT;
       mu_z2 = sqrt(1 - sinT2^2);
       mu_x2 = nWater * mu_x;
       mu_y2 = nWater * mu_y;
       delY2 = dAir * mu_y2/mu_z2;
       delX2 = dAir * mu_x2/mu_z2;
       
       ph_x = ph_x + delX1 + delX2;
       ph_y = ph_y + delY1 + delY2;       
       
       % find new position from air
       
       
       cos_rec_fov = cos(rec_fov(j)/2);
       
       distance = sqrt((ph_x-rx_x)^2 + (ph_y-rx_y)^2);                      % Euclidian distance to receiver center
       
       if ((distance <= radius) && (mu_z2 >= cos_rec_fov))   % Photon received
          power(j) = power(j) + exp(-a*total_rec_dist(i));          % Adjust power of received photons
          ph_cnt(j) = ph_cnt(j) + 1;                                % Increment number of received photons
          ph_power = exp(-a*total_rec_dist(i));
          
          % Record stats on incident angle
          angle_mean(j) = ((ph_cnt(j)-1)*angle_mean(j) + rec_loc_final(i,5)) / ph_cnt(j);   % E[theta] = (N*E[theta] + theta) / N+1 
          angle_sqravg = ((ph_cnt(j)-1)*angle_sqravg + rec_loc_final(i,5)^2)/ph_cnt(j);     % E[theta^2]
          % Record stats on distance traveled 
          power_mean(j) = ((ph_cnt(j)-1)*power_mean(j) + ph_power) / ph_cnt(j);
          power_sqravg = ((ph_cnt(j)-1)*power_sqravg +  ph_power^2) / ph_cnt(j);
          
          if ph_cnt(j) == 1
            angle_var(j) = (angle_sqravg - angle_mean(j)^2);    % Don't correct with N/N-1 since, it's N/0
            power_var(j) = (power_sqravg - power_mean(j)^2);
          else
            angle_var(j) = (ph_cnt(j)/(ph_cnt(j)-1))*(angle_sqravg - angle_mean(j)^2);    % Var[theta] = E[theta^2] - E[theta]^2  
            power_var(j) = (ph_cnt(j)/(ph_cnt(j)-1))*(power_sqravg - power_mean(j)^2);    % Var[P] = E[P^2] - E[P]^2  
          end
       end
       
   end
end

%power_mean.*ph_cnt - power