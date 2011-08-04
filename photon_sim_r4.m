% MC simulation for photon movementr
% V3 - changed propogation to z-axis, changed equations for updating
% position/angle
% V4 - Lots of changes to statistics calculations - mean and variance


% clear all
% clc

% Set the random stream seed to something ... wait for it ... random
RandStream.setDefaultStream ...
     (RandStream('mt19937ar','seed',sum(100*clock)));

if (isunix())
    userData = urlread('http://169.254.169.254/latest/user-data');
    if (~strcmp(userData,'autorun_sim'))
        %error('No autorun');
    else
        autorun = 1;
        cd('/home/wccox/Dropbox/WCC Research/mc');
        
        % Define these variables appropriately:
        mail = 'wccoxresearch@gmail.com'; %Your GMail email address
        load password

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

        sendmail('gallamine@gmail.com','Simulation started on AWS.')
    end
end

if (isunix()) 
    dataDir = '/home/wccox/';
    simDir = '/home/wccox/Dropbox/WCC Research/mc';
else
    dataDir = 'C:\Users\wccox\Documents\ThesisData\TankSimulations\RoundTwo';
    simDir = 'C:\Users\wccox\Dropbox\WCC Research\mc';
end

useVCL = 'false';                   % save output file to k drive
sendEmail = 'true';                % Send email at start and stop of simulation
saveOutput = 'true';                % Save the output data to a folder
ftpData = 'false';                  % FTP data back to FTP server at conclusion of simulation

num_photons = 1e6;                  % number of photons simulated per batch/group
num_sims = 500;                     % number of groups to simulate
n_water = 1.33;                     % index of refraction of water
n_window = 1.585;                   % index of refraction of polycarbonate

diverg = 0;

% [cdf_scatter,angle] = generate_scatter('measured','maalox_alan_orig');
% [cdf_scatter,angle] = generate_scatter('measured','petzold_maalox');
[cdf_scatter,angle] = generate_scatter('measured','widemann_maalox');
% [cdf_scatter,angle] = generate_scatter('measured','petzold_avg');
% [cdf_scatter,angle] = generate_scatter('measured','mie_1_micron');


% albedo = (c-a)/c;   % Water albedo is scattering coef./atten. coef. (b/c unitless)
albedo = 0.95;          % Albedo of Maalox (ranges from 0.8 to 0.95) - IF YOU CHANGE THIS, BE SURE TO CHANGE THE MINIMUM POWER VALUE!!!

c = 6.830601;
b = c * albedo;
a = c - b;

beamDiverg = 0.0015;                       %(0.01)*pi/180;           %degtorad(0.01);
beamWidth = 0.001;                  % 1.6 mm (half width). Hecht pg. 595

% beamDiverg = 0;
% beamWidth = 0;

rec_pos = [0,0];

sizeRecPos = size(rec_pos);
num_rx = sizeRecPos(1);

rec_aperture = ones(num_rx,1).*0.8;
rec_fov = ones(num_rx,1).* pi./2;

% rec_aperture = [0.0079, 0.0508, 0.0508]; % rec_aperture = ones(num_rx,1).*0.0508;  % 0.0508 m = 2 inches
% rec_fov = ones(num_rx,1).*0.314159;     % 0.314159 = 18 deg FOV; ;  % 0.0508 m = 2 inches
% rec_fov = [2.27, 0.0785, 0.1745];     % 0.314159 = 18 deg FOV, 2.27 = 130 deg, 0.0785 = 5 deg

receiver_z = 3.66;                     % Z position of the receiver (in meters)

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

attenLen = round(receiver_z*c);
                                        
scattering_events = ceil((c-a)*receiver_z*7)     %five times the scattering attenuation length
if scattering_events < 10
      scattering_events = 10;
end

photonDist = 0;         % Array to hold receiver distances of photons from receiver 
photonAngles = 0;
photonWeights = 0;

init_angle = 0;       % Point transmitter at receiver
init_angle2 = 0;       % Point transmitter at receiver


if (strcmp(saveOutput,'true'))
    foldername = sprintf('outputData-%s',datestr(now,'HH-MM-SS_yyyy-mm-dd'));
    mkdir(sprintf('%s/%s',dataDir,foldername));
    fid = fopen(sprintf('%s/%s/whatsHere.txt',dataDir,foldername),'w');
    fprintf(fid,'Experiment run on %s \n\r Atten. Length: %d \n\r Rx. dist: %d c: %d a: %d b: %d albedo: %d\n\r Num. photons: %d Num. sims %d',...
        foldername,attenLen,receiver_z,c,a,b,albedo,num_photons,num_sims);
    save(sprintf('%s/%s/simVariables.mat',dataDir,foldername));
    fclose(fid);
end

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


tStart = tic;

poolSize = 8;
if (matlabpool('size')==0)
    matlabpool('open','local',poolSize) 
end

totalPhotonsAtRxPlane = 0;
run_counter = 0;

recPosX = 0;
recPosY = 0;
finalPhotonPos = zeros(1,5);
finalPhotonDist = 0;
finalPhotonWeight = 0;

parfor simcount = 1:num_sims
    %%
    simcount
    [total_time(simcount), ...
    total_rec_power, ...
    total_rec_packets, ...
    rec_loc_final, ...
    total_rec_dist, ...
    rec_weights] = ...
    ...
    mc_func_r6(num_photons,...
    scattering_events,...
    c,...
    a,...
    receiver_z,...
    cdf_scatter,...
    angle,...
    init_angle,...
    init_angle2,...
    beamDiverg,beamWidth);

    totalPhotonsAtRxPlane = totalPhotonsAtRxPlane + total_rec_packets;
    
%     finalPhotonPos = [finalPhotonPos;rec_loc_final];
%     finalPhotonDist = [finalPhotonDist;total_rec_dist];
%     finalPhotonWeight = [finalPhotonWeight;rec_weights];
    
%======================= CODE FOR RECEIVER ================================    
%     [power,ph_cnt,angle_mean,angle_var,dist_mean,dist_var,weight_mean,weight_var,reflected,distances,angles,weights] = mc_rec_r4(a,rec_loc_final,total_rec_dist,rec_weights,rec_pos,rec_aperture,rec_fov,num_photons); 
%     
%     total_power = total_power + power';                             % Vectorized sum of the weights of received photons (sum received photons weights over all groupings)
%     total_photons = total_photons + ph_cnt';                        % Vectorized sum of number of photons
%     photonCount(simcount,:) = ph_cnt;                       % Store the total received photons/detector - used to weight statistics at the end
%      
%     photonDist = [photonDist distances];
%     photonAngles = [photonAngles angles];
%     photonWeights = [photonWeights weights];
% 
%     % Calculate MEAN and VARIANCE of the ANGLE 
%    
%     angleMean(simcount,:) = angle_mean;                        % Receiver angle means for each sub-simulation
%     angleVarSum = angleVarSum + (ph_cnt' - 1).*angle_var';
%     
%     % Calculate MEAN and VARIANCE of the DISTANCE 
% 
%     distMean(simcount,:) = dist_mean;
%     distVarSum = distVarSum + (ph_cnt' - 1).*dist_var';
%      
%     % Calculate MEAN and VARIANCE of the WEIGHT (each sub-simulation has a
%     % constant number of samples, so that's why the equation is different)
% 
%     weightMean(simcount,:) = weight_mean;
%     weightVarSum = weightVarSum + weight_var';
%     
%     reflec = reflected/total_rec_packets;
        
    if (strcmp(saveOutput,'true'))
        parsave(foldername,simcount,dataDir,rec_loc_final,total_rec_dist,rec_weights);
    end
end

% Check to make sure total_photons > 1, otherwise these values below will be NaN

% ======== CALCULATE SOME STATISTICS FROM OUTPUT OF RECEIVER ============== 
% totalMeanAngle = sum(photonCount.*angleMean,1)./total_photons';
% totalVarAngle = (1./(total_photons'-1)).*(angleVarSum' + sum(photonCount.*(angleMean - repmat(totalMeanAngle,size(angleMean,1),1)).^2,1));
% 
% totalMeanDist = sum(photonCount.*distMean,1)./total_photons';
% totalVarDist = (1./(total_photons'-1)).*(distVarSum' + sum(photonCount.*(distMean - repmat(totalMeanDist,size(distMean,1),1)).^2,1));
% 
% totalMeanWeight = sum(weightMean,1)./num_sims;
% totalVarWeight = ((num_photons-1)/(num_photons*num_sims - 1)).*weightVarSum' + ((num_photons)/(num_photons*num_sims - 1)).*sum(weightMean - repmat(totalMeanWeight,size(weightMean,1),1),1).^2;
% 


% disp('total_power/(num_photons*num_sims) = ');
% normRxPower = total_power./(num_photons*num_sims)               % Total received power / total transmitted photon packets SHOULD EQUAL total_mean_weight
% sd_dist = sqrt(totalVarDist);     

sim_time = toc(tStart) / 60

% ======== PLOT VARIOUS FIGURES (histograms of Brandon's data) ===========
% % Distance 
% binMin = 0;
% binMax = 0.40;
% numBins = 200;
% figure(4);
% hold on;
% 
% edges = [binMin:binMax/numBins:binMax];
% cnt = weightedhistc(photonDist,photonWeights,edges,'right');
% cnt = cnt(1:end-1);
% binDelta = ones(1,length(cnt)).*binMax/numBins;
% cntNorm = cnt./(pi.*(2.*edges(2:end).*binDelta - binDelta.^2));     % Remove bias from the annular ring area (area -> distance)
% cntNorm = cntNorm ./ cntNorm(1);                            % Normalize to peak value
% semilogy(edges(1:end-1),cntNorm,'.');
% 
% 
% % normalize based on the end value
% figure(1); hold on;
% edges = [binMin:binMax/numBins:binMax];
% cnt = weightedhistc(photonDist,photonWeights,edges,'right');
% cnt = cnt(1:end-1);
% binDelta = ones(1,length(cnt)).*binMax/numBins;
% cntNorm = cnt./(pi.*(2.*edges(2:end).*binDelta - binDelta.^2));     % Remove bias from the annular ring area (area -> distance)
% cntNorm = cntNorm ./ cntNorm(end);                            % Normalize to peak value
% semilogy(edges(1:end-1),cntNorm,'.');


if (matlabpool('size')>0)
    %matlabpool close
end
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

filename = sprintf('outputData-%s.mat',datestr(now,'HH-MM-SS_yyyy-mm-dd'));

if (strcmp(ftpData,'true'))
    % zip output data folder
    cd('~/');
    tar(foldername,foldername,'/home/wccox/')
    filename = sprintf('%s.tar',foldername);
    cd('/home/wccox/Dropbox/WCC Research/mc');
    dir = '~/';
    saveDataFTP(filename,dir); 
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
    %save output.mat
   
       
    % Define these variables appropriately:
    mail = 'wccoxresearch@gmail.com'; %Your GMail email address
    load password

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