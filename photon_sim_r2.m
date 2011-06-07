
clear all
clc

useVCL = 'false';
sendEmail = 'true';
saveOutput = 'true';

num_photons = 4e6;
num_sims = 140; 
n_water = 1.33;                     % index of refraction of water

%diverg = 0.01*pi/180;
diverg = 0;

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
% % Clear water maalox
% c = 0.1514;
% a = 0.03028;
% % Coastal water maalox
% c = 0.399;
% a = 0.0798;
% % Harbor water maalox
% c = 2.19;
% a = 0.438;

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

% c = 0.82;
% a = 0.164;
% c = 2.459;
% a = 0.4918;
% c = 3.552;
% a = 0.710;

[cdf_scatter,angle] = generate_scatter('measured','maalox_alan');


% albedo = (c-a)/c;   % Water albedo is scattering coef./atten. coef. (b/c unitless)
albedo = 0.8;          % Albedo of Maalox (ranges from 0.8 to 0.95)

c = 5.1;
b = c * albedo;
a = c - b;

% rec_pos = [0,0;
%            0.01,0;
%            0.03,0;
%            0.05,0;
%            0.07,0;
%            0.1,0;
%            0.15,0;
%            0.2,0;
%            0.25,0;
%            0.3,0];

rec_pos = [0,0; 0,0; 0,0];   % receivers at the same position
       
num_rx = length(rec_pos);
% rec_aperture = ones(num_rx,1).*0.0508;  % 0.0508 m = 2 inches
% rec_fov = ones(num_rx,1).*0.314159;     % 0.314159 = 18 deg FOV

rec_aperture = [0.0079, 0.0508, 0.0508]; % rec_aperture = ones(num_rx,1).*0.0508;  % 0.0508 m = 2 inches
% rec_fov = ones(num_rx,1).*0.314159;     % 0.314159 = 18 deg FOV; ;  % 0.0508 m = 2 inches
rec_fov = [2.27, 0.0785, 0.1745];     % 0.314159 = 18 deg FOV, 2.27 = 130 deg, 0.0785 = 5 deg


receiver_x = 3.66;                     % X position of the receiver (in meters)

% Diameter of aperture (in meters): 
% 0.0508m = 2"
% 0.0079m = OP 2 VIS diameter (7.9 mm)

% Field of view (in radians)
% 3 deg. -> 0.0523
% 4.5 deg -> 0.0785 rad
% 10 deg -> 0.1745
% 18 deg -> 0.314259 rad
% 95 deg -> 1.6581
% 130 deg -> 2.269 rad

                                        
scattering_events = ceil((c-a)*receiver_x*5)     %five times the scattering attenuation length
if scattering_events < 10
      scattering_events = 10;
end


init_angle = 0;       % Point transmitter at receiver
init_angle2 = 0;       % Point transmitter at receiver

total_time = zeros(num_sims,1);
total_rec_power = zeros(num_sims,1);
total_rec_packets = zeros(num_sims,1);
received_location = [];
travel_distance = [];

total_power = zeros(num_rx,1);
total_photons = zeros(num_rx,1);
total_variance_angle = zeros(num_rx,1);
total_mean_angle = zeros(num_rx,1);
total_variance_dist = zeros(num_rx,1);
total_mean_dist = zeros(num_rx,1);
total_mean_angle_sum = zeros(num_rx,1);
total_mean_dist_sum = zeros(num_rx,1);


% receiver_x = (1/c).*scattering_length;
% c = [13.66/num_sims:13.66/num_sims:13.66];
% c = [(3 .* ones(1,40)) (3.27 .* ones(1,50))]

% c = scattering_length./receiver_x;
% a = c.*(1-albedo);                      % albedo = b/c. a = c-b. a = c-c*albedo = c(1-albedo)

tStart = tic;

if (matlabpool('size')==0)
    matlabpool open local 8
end

totalPhotonsAtRxPlane = 0;
run_counter = 0;

parfor simcount = 1:num_sims
    %%
    simcount
    %%disp('Starting simulation: %d',simcount)
    [total_time, ...
    total_rec_power, ...
    total_rec_packets, ...
    rec_loc_final, ...
    total_rec_dist] = ...
    mc_func_r4(num_photons,...
    scattering_events,...
    c,...
    a,...
    receiver_x,...
    cdf_scatter,...
    angle,...
    init_angle,...
    init_angle2,...
    diverg);

    totalPhotonsAtRxPlane = totalPhotonsAtRxPlane + total_rec_packets;
    
    [power,ph_cnt,angle_mean,angle_var,distance_mean,distance_var] = mc_rec_r2(a,rec_loc_final,total_rec_dist,rec_pos,rec_aperture,rec_fov); 
    
    total_power = total_power + power';
    total_photons = total_photons + ph_cnt';
    total_variance_angle = total_variance_angle + angle_var';
    total_mean_angle_sum = total_mean_angle_sum + angle_mean';
    total_variance_dist = total_variance_dist + distance_var';
    total_mean_dist_sum = total_mean_dist_sum + distance_mean';
    run_counter = run_counter+1;

    
    parsave(power,ph_cnt,angle_mean,angle_var,distance_mean,distance_var,num_photons,simcount);

end

total_mean_angle = total_mean_angle_sum / num_sims;
total_mean_dist = total_mean_dist_sum / num_sims;

sim_time = toc(tStart) / 60


%     figure(1)
%     scatter(photon(:,1),photon(:,2),50,log10(photon(:,6)),'.')

if (matlabpool('size')>0)
    %matlabpool close
end

minutes_to_run = total_time/60

% %% Separate multiple simulation runs into one experiment
% k=1;
% sim_power(k) = total_rec_power(1);
% sim_count(k)=1;
% 
% for i=2:num_sims
%    if receiver_y(i) == receiver_y(i-1)
%         sim_power(k) = sim_power(k)+total_rec_power(i);
%         sim_count(k) = sim_count(k)+1;
%    else
%        k=k+1;
%        sim_power(k) = total_rec_power(i);
%        sim_count(k)=1;
%    end
% end
% sim_power_norm = sim_power./(sim_count.*num_photons)
% %%
% 
% i = 1;
% k = 1;
% total_packets_sim = zeros(1,length(sim_count));
% while i <= length(sim_count)
%     for j = 1:sim_count(i)
%         total_packets_sim(i) = total_packets_sim(i) + total_rec_packets(k);
%         k = k +1;
%     end
%     i = i+1;
% end
% 
% %%


% figure(8)                   % plot 3D histogram of photons on RX plane
% hist3(received_location)
% set(gcf,'renderer','opengl');
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');

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
sprintf('C = %d (1/m), A = %d (1/m).',c,a)
% sprintf('Receiver at %d, %d, %d (meters)',receiver_x, receiver_y, receiver_z)
% sprintf('Travel distance delta %d (m). Time of arrival delta %d (sec)', distance_delta, time_delta)
% sprintf('Time delta between histogram bins: %d (sec), %d (Hz)',T,bandwidth)

beep
beep
beep
beep

if (strcmp(useVCL,'true'))
    cd('K:\MC_data')
    save output.mat
end

if (strcmp(sendEmail,'true'))
    
    [ret, name] = system('hostname');   

    if ret ~= 0,
       if ispc
          name = getenv('COMPUTERNAME');
       else      
          name = getenv('HOSTNAME');      
       end
    end
    name = lower(name);
    
    subject = sprintf('Simulation on %s with %d photons, %d scattering events completed.', name, num_photons*num_sims, scattering_events)
    body = sprintf('C = %d (1/m), A = %d (1/m). Simulation took %d minutes to run.',c,a,minutes_to_run);
    save output.mat
   
       
    % Define these variables appropriately:
    mail = 'wccoxresearch@gmail.com'; %Your GMail email address
    password = ''; %Your GMail password

    % Then this code will set up the preferences properly:
    setpref('Internet','E_mail',mail);
    setpref('Internet','SMTP_Server','smtp.gmail.com');
    setpref('Internet','SMTP_Username',mail);
    setpref('Internet','SMTP_Password',password);
    props = java.lang.System.getProperties;
    props.setProperty('mail.smtp.auth','true');
    props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
    props.setProperty('mail.smtp.socketFactory.port','465');

    % Send the email. Note that the first input is the address you are sending the email to
    
    sendmail('gallamine@gmail.com','Simulation complete',[subject body],{'output.mat'}) 
end