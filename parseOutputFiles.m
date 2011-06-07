% Load output files saved from multiple simulations and combine the
% results. THIS NEEDS TO BE CHECKED!

simNumber = 10;
num_rx = 3;


total_power = zeros(num_rx,1);
total_photons = zeros(num_rx,1);
total_variance_angle = zeros(num_rx,1);
total_mean_angle = zeros(num_rx,1);
total_variance_dist = zeros(num_rx,1);
total_mean_dist = zeros(num_rx,1);
total_mean_angle_sum = zeros(num_rx,1);
total_mean_dist_sum = zeros(num_rx,1);

% Each output file contains:
% power,ph_cnt,angle_mean,angle_var,distance_mean,distance_var,num_photons,run_counter

for i = 1:simNumber
    load(sprintf('output%d.mat', i))
    
    total_power = total_power + power';
    total_photons = total_photons + ph_cnt';
    total_variance_angle = total_variance_angle + angle_var';
    total_mean_angle_sum = total_mean_angle_sum + angle_mean';
    total_variance_dist = total_variance_dist + distance_var';
    total_mean_dist_sum = total_mean_dist_sum + distance_mean';
    
    

end

total_mean_angle = total_mean_angle_sum / simNumber;
total_mean_dist = total_mean_dist_sum / simNumber;
totalTxPhotons = num_photons*simNumber;