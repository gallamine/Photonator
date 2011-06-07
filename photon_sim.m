
clear all
clc


num_photons = 2e6;
scattering_events = 15; 
n_water = 1.33;                     % index of refraction of water

% % harbor water
% c = 2.190;
% a = 0.366;
% [cdf_scatter,angle] = generate_scatter('calc','harbor');
% % coastal water
% c = 0.22 + 0.179;
% a = 0.179;
% [cdf_scatter,angle] = generate_scatter('calc','coastal');
% % Clear water
% c = 0.0374 + 0.114;
% a = 0.114;
% [cdf_scatter,angle] = generate_scatter('calc','clear');
% Maalox
% c = 0.275;
% a = 0.055;
% c = 0.5;
% a = 0.1;
% c = 1.0;
% a = 0.2;
% c = 1.475;
% a = 0.295;
% c = 1.92;
% a = 0.384;
% c = 2.25;
% a = 0.45;
% c = 2.675;
% a = 0.535;
% c = 3.2;
% a = 0.64;
[cdf_scatter,angle] = generate_scatter('measured','maalox_alan');

% scattering_length = zeros(45,1);
% scattering_length(1:3:end) = [1:2:30];
% scattering_length = conv([1 1 1],scattering_length);

% scattering_length = [1;
%     3;
%     5;5;
%     7;7;
%     9;9;9;
%     11;11;11;11;
%     13;13;13;13;13;13;13;
%     15;15;15;15;15;15;15;
%     17;17;17;17;17;17;17;17;17;17;
%     19;19;19;19;19;19;19;19;19;19;19;19;19];

% albedo = (c-a)/c;   % Water albedo is scattering coef./atten. coef. (b/c unitless)
albedo = 0.95;          % Albedo of Maalox

c = 1;
b = c * albedo;
a = c - b;

h = waitbar(0.0,'Please wait...','CreateCancelBtn','stop=true; delete(h); clear h');  
set(h,'Name','Simulation Progress');  

receiver_x = 3.66;                     % X position of the receiver (in meters)
% receiver_y = [0;
%              0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01;
%              0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;0.05;
% 
%              0.3;0.3;0.3;0.3;0.3;0.3;0.3;0.3;0.3;0.3;0.3; 0.3;0.3;
%              0.3;0.3;0.3;0.3;0.3;0.3;0.3;0.3;0.3;0.3;0.3; 0.3;0.3;     
%              0.40;0.40;0.40;0.40;0.40;0.40;0.40;0.40;0.40;0.40;0.40;
%              0.40;0.40;0.40;0.40;0.40;0.40;0.40;0.40;0.40;0.40;0.40;0.40;0.40;0.40;0.40;];                     % Y position of the receiver (in meters)
receiver_y = ones(8,1)*0.00;
receiver_z = 0;                     % Z position of the receiver(in meters)
aperture = 0.0079;                       % Diameter of aperture (in meters)

fov = 2.269;                           % Field of view (in radians)
                                        % 3 deg. -> 0.0523
                                        % 10 deg -> 0.1745
                                        % 95 deg -> 1.6581
                                        % 130 deg -> 2.269 rad

% init_angle = atan(receiver_y/receiver_x)       % Point transmitter at receiver
% init_angle2 = atan(receiver_z/receiver_x)       % Point transmitter at receiver
init_angle = 0;       % Point transmitter at receiver
init_angle2 = 0;       % Point transmitter at receiver

num_sims = length(receiver_y);                                   % How many times to run the simulations
total_time = zeros(num_sims,1);
total_rec_power = zeros(num_sims,1);
total_rec_packets = zeros(num_sims,1);
received_location = [];
travel_distance = [];

% receiver_x = (1/c).*scattering_length;
% c = [13.66/num_sims:13.66/num_sims:13.66];
% c = [(3 .* ones(1,40)) (3.27 .* ones(1,50))]

% c = scattering_length./receiver_x;
% a = c.*(1-albedo);                      % albedo = b/c. a = c-b. a = c-c*albedo = c(1-albedo)

tStart = tic;

if (matlabpool('size')==0)
    matlabpool open local 8
end

parfor simcount = 1:num_sims
    %waitbar(simcount/num_sims,h,['Simulation ' num2str(simcount)]); 
    [total_time(simcount),total_rec_power(simcount),total_rec_packets(simcount),rec_loc_final,total_rec_dist] = ...
    mc_func_r3(num_photons,scattering_events,c,a,receiver_x,receiver_y(simcount),receiver_z,...
    aperture,fov,cdf_scatter,angle,init_angle,init_angle2);

    received_location = vertcat(received_location,rec_loc_final);
    travel_distance = vertcat(travel_distance,total_rec_dist);
    
end

sim_time = toc(tStart) / 60
close(h)

%     figure(1)
%     scatter(photon(:,1),photon(:,2),50,log10(photon(:,6)),'.')

if (matlabpool('size')>0)
    %matlabpool close
end

minutes_to_run = total_time/60

%% Separate multiple simulation runs into one experiment
k=1;
sim_power(k) = total_rec_power(1);
sim_count(k)=1;

for i=2:num_sims
   if receiver_y(i) == receiver_y(i-1)
        sim_power(k) = sim_power(k)+total_rec_power(i);
        sim_count(k) = sim_count(k)+1;
   else
       k=k+1;
       sim_power(k) = total_rec_power(i);
       sim_count(k)=1;
   end
end
sim_power_norm = sim_power./(sim_count.*num_photons)
%%

i = 1;
k = 1;
total_packets_sim = zeros(1,length(sim_count));
while i <= length(sim_count)
    for j = 1:sim_count(i)
        total_packets_sim(i) = total_packets_sim(i) + total_rec_packets(k);
        k = k +1;
    end
    i = i+1;
end

%%


figure(8)                   % plot 3D histogram of photons on RX plane
hist3(received_location)
set(gcf,'renderer','opengl');
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');

% figure(10)                              % Plot histogram of time-of-arrival vs. power
% [N,X] = hist(travel_distance,100);
% pwr_rx_hist = N.*exp(-a.*X);
% stem(X/(3e8/n_water),pwr_rx_hist/sum(total_rec_power))   
% xlabel('Time of arrival (sec)')
% ylabel('Percentage of power')

% distance_delta = max(X) - min(X);
% time_delta = distance_delta/(3e8/n_water);
% T = mean(X(2:end) - X(1:end-1))/(3e8/n_water); % Effective "sampling rate" of the histogram
% bandwidth = 1/T;                                % Normalized frequency in Hz
% 
%  figure(7)
% freqz(N/sum(total_rec_packets),[1],512,bandwidth)      %% plot a frequency response from the histogram data

% figure(9)                           %% Scatter plot of photons on RX plane
% scatter(rec_loc(:,1),rec_loc(:,2))

% figure(1)
% scatter(photon(:,1),photon(:,2),50,log10(photon(:,6)),'.')
% xlim([-20 150]);
% ylim([-120 120]);
% xlabel('x-axis (m)')
% xlabel('y-axis (m)')
% 
% figure(4)
% scatter3(photon(:,1),photon(:,2),photon(:,3),50,log10(photon(:,6)),'.')

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

sprintf('Simulation on DATE with %d photons, %d scattering events.', num_photons*num_sims, scattering_events)
sprintf('C = %d (1/m), A = %d (1/m). FOV = %d (radians). Aperture = %d (m).',c,a,fov,aperture)
sprintf('Receiver at %d, %d, %d (meters)',receiver_x, receiver_y, receiver_z)
% sprintf('Travel distance delta %d (m). Time of arrival delta %d (sec)', distance_delta, time_delta)
% sprintf('Time delta between histogram bins: %d (sec), %d (Hz)',T,bandwidth)

beep
beep
beep
beep