% standaloneReceiver
% This program reads in stored datafiles of photons on the receiver plane,
% and passes the data to the receiver function for processing. This allows
% the photon movement (to the Rx plane) and the "reception" to be
% separated.

if (isunix()) 
    dataDirectory = '/home/wccox/outputData-21-25-32_2011-08-11';
else
    disp('Windows');
    dataDirectory = 'C:\Users\wccox\Documents\ThesisData\TankSimulations\RoundTwo\outputData-07-47-36_2011-08-05'; 
end

% a = 0.396;
% num_photons = 1e6;
load(sprintf('%s/simVariables.mat',dataDirectory));

file_list = dir(fullfile(dataDirectory,'*outputData*.mat'));

num_sims = size(file_list,1);

% % rec_fov = [3;6;18;27;45;90;130;180;3;6;18;27;45;90;130;180].*pi./180;
% rec_fov = [1;4;5;6;1;4;5;6].*pi./180;
% sizeRecPos = size(rec_fov);
% num_rx = sizeRecPos(1);
% 
% rec_pos = zeros(num_rx,2);
% % rec_aperture = [ones(num_rx/2,1).*0.0508; ones(num_rx/2,1).*0.1016];
% rec_aperture = [ones(num_rx/2,1).*0.0508; ones(num_rx/2,1).*0.0254];


rec_fov = 26*pi/180;                
num_rx = 1;
rec_pos = [0,0];
rec_aperture = 0.0508;


angleVarSum = zeros(num_rx,1);
angleMean = zeros(num_sims,num_rx);
distVarSum = zeros(num_rx,1);
distMean = zeros(num_sims,num_rx);
weightVarSum = zeros(num_rx,1);
weightMean = zeros(num_sims,num_rx);

total_power = zeros(num_rx,1);
total_photons = zeros(num_rx,1);
total_rec_packets = zeros(num_sims,1);


tic;
parfor simcount = 1:num_sims
    %%
    simcount
    S = load(sprintf('%s/%s',dataDirectory,file_list(simcount).name));
%     foldername = S.varargin{1}
%     simcount =  S.varargin{2}
%     dataDirectory = S.varargin{3}           <- this changed after 1st data set!
%     rec_loc_final = S.varargin{4}
%     total_rec_dist = S.varargin{5}
%     rec_weights = S.varargin{6}

    
%======================= CODE FOR RECEIVER ================================    
    [power,ph_cnt,angle_mean,angle_var,dist_mean,dist_var,weight_mean,weight_var,reflected] ...
        = mc_rec_r5(a,S.varargin{4},S.varargin{5},S.varargin{6},rec_pos,rec_aperture,rec_fov,num_photons); 
    
    total_power = total_power + power';                             % Vectorized sum of the weights of received photons (sum received photons weights over all groupings)
    total_photons = total_photons + ph_cnt';                        % Vectorized sum of number of photons
    photonCount(simcount,:) = ph_cnt;                       % Store the total received photons/detector - used to weight statistics at the end
     
%     photonDist = [photonDist distances];
%     photonAngles = [photonAngles angles];
%     photonWeights = [photonWeights weights];

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
    
    %reflec = reflected/total_rec_packets;

end

time = toc/60


% ======== CALCULATE SOME STATISTICS FROM OUTPUT OF RECEIVER ============== 
totalMeanAngle = sum(photonCount.*angleMean,1)./total_photons';
totalVarAngle = (1./(total_photons'-1)).*(angleVarSum' + sum(photonCount.*(angleMean - repmat(totalMeanAngle,size(angleMean,1),1)).^2,1));

totalMeanDist = sum(photonCount.*distMean,1)./total_photons';
totalVarDist = (1./(total_photons'-1)).*(distVarSum' + sum(photonCount.*(distMean - repmat(totalMeanDist,size(distMean,1),1)).^2,1));

totalMeanWeight = sum(weightMean,1)./num_sims;
totalVarWeight = ((num_photons-1)/(num_photons*num_sims - 1)).*weightVarSum' + ((num_photons)/(num_photons*num_sims - 1)).*sum(weightMean - repmat(totalMeanWeight,size(weightMean,1),1),1).^2;

varWeightPop = totalVarWeight./total_photons';

disp('total_power/(num_photons*num_sims) = ');
normRxPower = total_power./(num_photons*num_sims)               % Total received power / total transmitted photon packets SHOULD EQUAL total_mean_weight
sd_dist = sqrt(totalVarDist);  
