% VERSION LOG
% Rev 1 - first creation, record power and count for each receiver
% Rev 2 - add statistics, mean and variance for power, angle(s)

function [power,ph_cnt,angle_mean,angle_var,distance_mean,distance_var] = mc_rec_r2(a,rec_loc_final,total_rec_dist,rec_pos,rec_aperture,rec_fov)

num_rx = length(rec_pos);   % total number of receivers
num_photons = length(rec_loc_final);    % total number of received photons

ph_cnt = zeros(1,num_rx);
power = zeros(1,num_rx);

angle_mean = zeros(1,num_rx);   % mean of the incident angle for each Rx
angle_var = zeros(1,num_rx);    % variance of incident angle for each Rx
distance_mean = zeros(1,num_rx);    % Mean distance traveled for each Rx
distance_var = zeros(1,num_rx);     % Variance of dist traveled for each Rx

for j = 1:num_rx
   rx_y = rec_pos(j,1);
   rx_z = rec_pos(j,2);
   radius = rec_aperture(j)/2;
   angle_sqravg = 0;
   dist_sqravg = 0;
   
   %cos_rec_fov = cos(rec_fov(j)/2);
   
   for i = 1:num_photons
       ph_y = rec_loc_final(i,1);
       ph_z = rec_loc_final(i,2);
       distance = sqrt((ph_y-rx_y)^2 + (ph_z-rx_z)^2);                      % Euclidian distance to receiver center
       
       if ((distance <= radius) && (abs(rec_loc_final(i,3)) <= rec_fov(j)/2))   % Photon received
          power(j) = power(j) + exp(-a*total_rec_dist(i));          % Adjust power of received photons
          ph_cnt(j) = ph_cnt(j) + 1;                                % Increment number of received photons
          
          % Record stats on incident angle
          angle_mean(j) = ((ph_cnt(j)-1)*angle_mean(j) + rec_loc_final(i,3)) / ph_cnt(j);   % E[theta] = (N*E[theta] + theta) / N+1 
          angle_sqravg = ((ph_cnt(j)-1)*angle_sqravg + rec_loc_final(i,3)^2)/ph_cnt(j);     % E[theta^2]
          angle_var(j) = (ph_cnt(j)/(ph_cnt(j)-1))*(angle_sqravg - angle_mean(j)^2);    % Var[theta] = E[theta^2] - E[theta]^2
          % Record stats on distance traveled 
          distance_mean(j) = ((ph_cnt(j)-1)*distance_mean(j) + total_rec_dist(i)) / ph_cnt(j);
          dist_sqravg = ((ph_cnt(j)-1)*dist_sqravg +  total_rec_dist(i)^2) / ph_cnt(j);
          distance_var(j) = (ph_cnt(j)/(ph_cnt(j)-1))*(dist_sqravg - distance_mean(j)^2);
          
       end
       
   end
end