% Function to save output data from inside of a 'parfor' loop. This is
% hacky, but required.

% function parsave(power,ph_cnt,angle_mean,angle_var,distance_mean,distance_var,weight_mean,weight_var,num_photons,run_counter,reflected)
%   
%     save(sprintf('output%d.mat', run_counter),'power',...
%       'ph_cnt',...
%       'angle_mean','angle_var','distance_mean','distance_var',...
%       'weight_mean','weight_var',...
%       'num_photons','reflected','run_counter');

function parsave(varargin)

    
    filename = sprintf('%04d-%s.mat',varargin{2},varargin{1});   % output file name, sim #, folder name.mat
    %cd(varargin{1});
    save(sprintf('%s/%s/%s',varargin{3},varargin{1},filename),'varargin');
    %cd ..

end
