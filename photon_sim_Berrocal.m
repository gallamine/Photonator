% MC simulation for photon movementr
% V3 - changed propogation to z-axis, changed equations for updating
% position/angle
% V4 - Lots of changes to statistics calculations - mean and variance

% clear all
% clc

useVCL = 'false';
sendEmail = 'true';
saveOutput = 'false';

num_photons = 1e6;
num_sims = 8; 
n_water = 1.33;                     % index of refraction of water
n_window = 1.585;                   % index of refraction of polycarbonate

diverg = 0;

[cdf_scatter,angle] = generate_scatter('measured','mie_1_micron');


% albedo = (c-a)/c;   % Water albedo is scattering coef./atten. coef. (b/c unitless)
albedo = 1;          % Albedo of Maalox (ranges from 0.8 to 0.95) - IF YOU CHANGE THIS, BE SURE TO CHANGE THE MINIMUM POWER VALUE!!!

c = 500;
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

rec_pos = [0,0];   % receivers at the same position (on the x/y plane)
       
sizeRecPos = size(rec_pos);
num_rx = sizeRecPos(1);   % total number of receivers
% rec_aperture = ones(num_rx,1).*0.0508;  % 0.0508 m = 2 inches
% rec_fov = ones(num_rx,1).*0.314159;     % 0.314159 = 18 deg FOV

rec_aperture = [0.01]; % rec_aperture = ones(num_rx,1).*0.0508;  % 0.0508 m = 2 inches
% rec_fov = ones(num_rx,1).*0.314159;     % 0.314159 = 18 deg FOV; ;  % 0.0508 m = 2 inches
rec_fov = [0.296];     % 0.314159 = 18 deg FOV, 2.27 = 130 deg, 0.0785 = 5 deg, 0.296 = 8.5*2


receiver_z = 0.0050000001;                     % Z position of the receiver (in meters)

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

                                        
scattering_events = ceil((c-a)*receiver_z*7)     %five times the scattering attenuation length
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


angleVarSum = zeros(num_rx,1);
angleMean = zeros(num_sims,num_rx);
distVarSum = zeros(num_rx,1);
distMean = zeros(num_sims,num_rx);
weightVarSum = zeros(num_rx,1);
weightMean = zeros(num_sims,num_rx);

allWeights = 0;
allAngles = 0;
allDist = 0;

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

recPosX = 0;
recPosY = 0;

parfor simcount = 1:num_sims
    %%
    simcount
    %%disp('Starting simulation: %d',simcount)
    [total_time(simcount), ...
    total_rec_power, ...
    total_rec_packets, ...
    rec_loc_final, ...
    total_rec_dist, ...
    rec_weights] = ...
    ...
    mc_func_Berrocal(num_photons,...
    scattering_events,...
    c,...
    a,...
    receiver_z,...
    cdf_scatter,...
    angle,...
    init_angle,...
    init_angle2,...
    diverg);

    totalPhotonsAtRxPlane = totalPhotonsAtRxPlane + total_rec_packets;
    
    [power,ph_cnt,angle_mean,angle_var,dist_mean,dist_var,weight_mean,weight_var,reflected,ph_x,ph_y] = mc_rec_Berrocal(a,rec_loc_final,total_rec_dist,rec_weights,rec_pos,rec_aperture,rec_fov,num_photons); 
    
    recPosX = [recPosX ph_x];
    recPosY = [recPosY ph_y];
    
    total_power = total_power + power';                             % Vectorized sum of the weights of received photons (sum received photons weights over all groupings)
    total_photons = total_photons + ph_cnt';                        % Vectorized sum of number of photons
    photonCount(simcount,:) = ph_cnt;                       % Store the total received photons/detector - used to weight statistics at the end
    
%     allWeights = [allWeights weights];
%     allAngles = [allAngles angles];
%     allDist = [allDist distances];
    
    
    % Calculate MEAN and VARIANCE of the ANGLE 
   
    angleMean(simcount,:) = angle_mean;                        % Receiver angle means for each sub-simulation
    angleVarSum = angleVarSum + (ph_cnt' - 1).*angle_var';
    
    % Calculate MEAN and VARIANCE of the DISTANCE 

    distMean(simcount,:) = dist_mean;
    distVarSum = distVarSum + (ph_cnt' - 1).*dist_var';
     
    % Calculate MEAN and VARIANCE of the WEIGHT (each sub-simulation has a
    % constant number of samples, so that's why the equation is different)

    weightMean(simcount,:) = weight_mean;
    weightVarSum = weightVarSum + weight_var';
    
    reflec = reflected/total_rec_packets;
        
    if (strcmp(saveOutput,'true'))
        parsave(power,ph_cnt,angle_mean,angle_var,dist_mean,dist_var,weight_mean,weight_var,num_photons,simcount,reflected);
    end
end

% Check to make sure total_photons > 1, otherwise these values below will be NaN

totalMeanAngle = sum(photonCount.*angleMean,1)./total_photons';
totalVarAngle = (1./(total_photons'-1)).*(angleVarSum' + sum(photonCount.*(angleMean - repmat(totalMeanAngle,size(angleMean,1),1)).^2,1));

totalMeanDist = sum(photonCount.*distMean,1)./total_photons';
totalVarDist = (1./(total_photons'-1)).*(distVarSum' + sum(photonCount.*(distMean - repmat(totalMeanDist,size(distMean,1),1)).^2,1));

totalMeanWeight = sum(weightMean,1)./num_sims;
totalVarWeight = ((num_photons-1)/(num_photons*num_sims - 1)).*weightVarSum' + ((num_photons)/(num_photons*num_sims - 1)).*sum(weightMean - repmat(totalMeanWeight,size(weightMean,1),1),1).^2;

% totalVarWeight2 = var([allWeights(2:end) zeros(1,(num_photons*num_sims-length(allWeights(2:end))))])                        % allWeights has a leading 0
% totalVarAngle2 = var(allAngles(2:end))
% totalVarDist2 = var(allDist(2:end))


disp('total_power/(num_photons*num_sims) = ');
normRxPower = total_power./(num_photons*num_sims)               % Total received power / total transmitted photon packets SHOULD EQUAL total_mean_weight
sd_dist = sqrt(totalVarDist);     

sim_time = toc(tStart) / 60



figure(4);
cnt = hist3([recPosY' recPosX'],'Edges',{-0.005:5e-5:0.005 -0.005:5e-5:0.005});
finalPdf = cnt./max(max(cnt));
imshow(finalPdf,[0 1])
colormap(jet);




if (matlabpool('size')>0)
    %matlabpool close
end

%minutes_to_run = total_time/60


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
    body = sprintf('C = %d (1/m), A = %d (1/m). Simulation took %d minutes to run.',c,a,sim_time);
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
    
    sendmail('gallamine@gmail.com','Simulation complete',[subject body]) 
end