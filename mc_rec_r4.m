% VERSION LOG
% Rev 1 - first creation, record power and count for each receiver
% Rev 2 - add statistics, mean and variance for power, angle(s)
% Rev 3 - changed receiver plane to x/y instead of y/z - photons move along
% the z-axis now
% Rev 3b/4 - 5/3/11 - added window affects, changed power calculations

function [power,ph_cnt,angle_mean,angle_var,dist_mean,dist_var,weight_mean,weight_var,reflected,distances,angles,weights] = mc_rec_r4(a,rec_loc_final,total_rec_dist,rec_weights,rec_pos,rec_aperture,rec_fov,numTxPhotons)

sizeRecPos = size(rec_pos);
num_rx = sizeRecPos(1);   % total number of receivers

num_photons = length(rec_loc_final);    % total number of received photons

ph_cnt = zeros(1,num_rx);
power = zeros(1,num_rx);

angle_mean = zeros(1,num_rx);   % mean of the incident angle for each Rx
angle_var = zeros(1,num_rx);    % variance of incident angle for each Rx
dist_mean = zeros(1,num_rx);    % Mean distance traveled for each Rx

dist_var = zeros(1,num_rx);     % Variance of dist traveled for each Rx
weight_mean = zeros(1,num_rx);    % Mean weight for each Rx
weight_var = zeros(1,num_rx);     % Variance of weight for each Rx

distances =  zeros(1,num_photons);       % Distances from photon to 0,0 point - Initailize to full size, crop at end
angles =  zeros(1,num_photons);          % Angle of received photon
weights =  zeros(1,num_photons);         % Weight of received photon

dWindow = 0.00635;              % 0.25 inches
dAir = 0.127;                   % 5 inches
nWater = 1.33;
nWindow = 1.585;                % polycarbonate
nAir = 1;
reflected = 0;                  % number of photons past critical angle of window
nWaterWindow = nWater/nWindow;
nWindowAir = nWindow/nAir;

critAngSine = nAir/nWater;             % sine T1 = nair/nwater sine T3

if (num_rx > 1)
    disp('You need to change how the distance array is stored!');
end
   
for i = 1:num_photons                   % iterate over every photon on receiver plane
   ph_x = rec_loc_final(i,1);
   ph_y = rec_loc_final(i,2);
   mu_x = rec_loc_final(i,3);
   mu_y = rec_loc_final(i,4); 
   mu_z = rec_loc_final(i,5); 
   
   % Find new position from window
   sinT = sqrt(1 - mu_z^2);
   
%    if sinT >= critAngSine           % if the photon's incident angle > critical angle, then skip this photon
%        reflected = reflected + 1;
%        continue;
%    end
%    
%    sinT1 =  nWaterWindow*sinT;          % 1.33/1.585 * sqrt(1 - uz^2)
%    mu_z1 = sqrt(1-sinT1^2);             % this can be simplified ... 
%    mu_x1 = nWaterWindow * mu_x;
%    mu_y1 = nWaterWindow * mu_y;
%    delY1 = dWindow * mu_y1/mu_z1;
%    delX1 = dWindow * mu_x1/mu_z1;
%     % find new position from air
%    sinT2 = nWater*sinT;                     % this can be simplified ... 
%    mu_z2 = sqrt(1 - sinT2^2);               
%    mu_x2 = nWater * mu_x;
%    mu_y2 = nWater * mu_y;
%    delY2 = dAir * mu_y2/mu_z2;
%    delX2 = dAir * mu_x2/mu_z2;
% 
%    ph_x = ph_x + delX1 + delX2;
%    ph_y = ph_y + delY1 + delY2;    
   
   for j = 1:num_rx                     % iterate over every receiver 
       rx_x = rec_pos(j,1);
       rx_y = rec_pos(j,2);
       radius = rec_aperture(j)/2;          % 1/2 diameter of receiver
       angle_sqravg = 0;
       dist_sqravg = 0;   
       
       cos_rec_fov = cos(rec_fov(j)/2);                             % cos(fov/2) to compare with photon's incident angle
       
       distance = sqrt((ph_x-rx_x)^2 + (ph_y-rx_y)^2);              % Euclidian distance to receiver center
       
       if ((distance <= radius) && (mu_z >= cos_rec_fov))          % Photon received
          power(j) = power(j) + rec_weights(i);                     % Adjust power of received photons
          ph_cnt(j) = ph_cnt(j) + 1;                                % Increment number of received photons
          ph_dist = total_rec_dist(i);
          
          angle_delta = rec_loc_final(i,5) - angle_mean(j);
          dist_delta = ph_dist - dist_mean(j);
          weight_delta = rec_weights(i) - weight_mean(j);
          
%           if j == 1
%               weights(ph_cnt(j)) = rec_weights(i);
%               angles(ph_cnt(j)) = rec_loc_final(i,5);
%               distances(ph_cnt(j)) = ph_dist;
%           end
          
          % Record stats on incident angle
          angle_mean(j) = angle_mean(j) + (angle_delta) / ph_cnt(j);   % E'[uz] = E[uz] + (uz - E[uz]) / N+1 - See Knuth's Art of Comp. Programming & http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
          % Record stats on distance traveled 
          dist_mean(j) = dist_mean(j) + (dist_delta) / ph_cnt(j);
          % Record stats on photon weight
          weight_mean(j) = weight_mean(j) + (weight_delta)/ph_cnt(j);       % E'[D] = E[D] + (D - E[D]) / N+1
           
    
          angle_var(j) =  angle_var(j) + angle_delta*(rec_loc_final(i,5) - angle_mean(j));    % Var[uz] = 1/N-1 *(Var[uz] + (uz - E[uz]_N-1)(uz - E[uz]_N))  
          dist_var(j) = dist_var(j) + dist_delta*(ph_dist - dist_mean(j));    
          weight_var(j) = weight_var(j) + weight_delta*(rec_weights(i) - weight_mean(j));
          
          if (num_rx == 1)
              distances(ph_cnt(j)) = distance;
              angles(ph_cnt(j)) = mu_z;
              weights(ph_cnt(j)) = rec_weights(i);
          end
       end       
   end
end



for j = 1:num_rx
    
    if (num_rx == 1)
        distances = distances(1:ph_cnt(j));
        angles = angles(1:ph_cnt(j));
        weights = weights(1:ph_cnt(j));
    end
    
    weight_var(j) = (weight_var(j) + weight_mean(j).^2*(ph_cnt(j).*(numTxPhotons - ph_cnt(j))/numTxPhotons)) / (numTxPhotons - 1);
    weight_mean(j) = ph_cnt(j)*weight_mean(j) / numTxPhotons;
    
    if weight_mean(j) ~= power(j)/numTxPhotons
        disp('Problem counting up the weight mean');
        weight_mean(j) - power(j)/numTxPhotons
        weight_mean(j)
        power(j)
    end
    
%     if ph_cnt(j) > 1
%         angle_var(j) =  (1./(ph_cnt(j)-1)).*angle_var(j);
%         dist_var(j) = (1./(ph_cnt(j)-1)).*dist_var(j);
%     end
%     
%     if isnan(angle_var(j))
%         disp('angle_var(j) is NaN (in correction loop)');
%     end
%     if isnan(dist_var(j))
%         disp('dist_var(j) is NaN (in correction loop)');
%     end
%     if isnan(weight_var(j))
%         disp('weight_var(j) is NaN (in correction loop)');
%     end
%     if isinf(angle_var(j))
%         disp('angle_var(j) is Inf');
%     end
%     if isinf(dist_var(j))
%         disp('dist_var(j) is Inf');
%     end
%     if isinf(weight_var(j))
%         disp('weight_var(j) is Inf');
%     end

end



%dist_mean.*ph_cnt - power