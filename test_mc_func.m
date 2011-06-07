%%


num_photons = 1e3;
scattering_events = 1; 
n_water = 1.33;                     % index of refraction of water

c = 2.675;
a = 0.535;

[cdf_scatter,angle] = generate_scatter('measured','maalox_alan');

receiver_x = 0.1;                     % X position of the receiver (in meters)

init_angle = 0; %atan(3.66/receiver_x)      
init_angle2 = 0     

rec_pos = [0,0;1,1];    % y,z location of receiver
rec_aperture = [3;3];   % aperture diamter in meters
rec_fov = [90;20].*(pi/180);    % field of view (in degrees)

diverg = 0.01*pi/180;

[total_time,total_rec_power,total_rec_packets,rec_loc_final,total_rec_dist] = ...
    mc_func_r4(num_photons,scattering_events,c,a,receiver_x,cdf_scatter,angle,init_angle,init_angle2,diverg);
%%
[power,ph_cnt] = mc_rec_r1(a,rec_loc_final,total_rec_dist,rec_pos,rec_aperture,rec_fov);      
%%
% scatter(rec_loc_final(:,1),rec_loc_final(:,2))
% hist3([rec_loc_final(:,1),rec_loc_final(:,2)],[20,20])

%%