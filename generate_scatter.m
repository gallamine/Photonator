function [cdf_scatter,petzold_angle_rad] = generate_scatter(func_type,water_cond,g)

if strcmp(func_type,'')
    func_type = 'measured';
end

if strcmp(water_cond,'')
    water_cond = 'clear';
end

%func_type = 'calc'      % 'measured' or 'calc'

% %     g = 0.924;
% % GENERAL USAGE:
% % [cdf_scatter,angle] = generate_scatter(0.924);
% %     angle1 = 0:0.006:pi;
% %     angle2 = -pi:0.006:0;
%     angle2 = -0.25:0.001:(2*pi-0.25);
% %     angle = -pi:0.006:pi;
% 
% 
% %     scatter1 = (1/(4*pi))*((1-g^2)./(1+g^2 - 2*g.*cos(angle1)).^1.5);
%     scatter2 = (1/(4*pi))*((1-g^2)./(1+g^2 - 2*g.*cos(abs(angle2))).^1.5);
% %     scatter = [scatter2 scatter1]; 
%     scatter = scatter2;
%     scatter_norm = scatter ./ sum(scatter);
%     cdf_scatter = cumsum(scatter_norm);
% 
% load petzold_data_orig
% petzold_angle = [0.1:0.05:180]';
% cdf_scatter = cumsum(scatter_water.*sind(angle_water));       
% cdf_scatter = interp1(angle_water,cdf_scatter,petzold_angle);
% cdf_scatter = cdf_scatter ./ max(cdf_scatter);
% petzold_angle_rad = petzold_angle.*pi./180;

if strcmp(func_type,'measured')
    
    if strcmp(water_cond,'petzold_avg')

        load petzold_data_orig
        %petzold_angle = [0.1:0.05:180]';

        cdf_scatter = cumtrapz(angle_water,scatter_water.*sind(angle_water));
        %cdf_scatter = interp1(angle_water,cdf_scatter,petzold_angle);
        cdf_scatter = cdf_scatter ./ max(cdf_scatter);
        petzold_angle_rad = angle_water.*pi./180;
        
    elseif strcmp(water_cond,'maalox_alan')
       
        load maalox_alan
        %load maalox_alan % these are diretly measured, and aren't extrapolated for small angles
        
%         petzold_angle = [maalox_angle(1:99); [4.30:0.05:180]'];   %only use if not extrapolating values
        %petzold_angle = [maalox_angle(1:144); [4.30:0.05:180]'];

        cdf_scatter = cumtrapz(maalox_angle,maalox_vsf.*sind(maalox_angle));       
        %cdf_scatter = interp1(maalox_angle,cdf_scatter,petzold_angle);
        cdf_scatter = cdf_scatter ./ max(cdf_scatter);
        %petzold_angle_rad = petzold_angle.*pi./180;
        petzold_angle_rad = maalox_angle.*pi./180;
        
    elseif strcmp(water_cond,'maalox_alan_orig')
               
        load maalox_alan_orig
        cdf_scatter = cumtrapz(maalox_angle,maalox_vsf.*sind(maalox_angle));       
        %cdf_scatter = interp1(maalox_angle,cdf_scatter,petzold_angle);
        cdf_scatter = cdf_scatter ./ max(cdf_scatter);
        %petzold_angle_rad = petzold_angle.*pi./180;
        petzold_angle_rad = maalox_angle.*pi./180;
        
    elseif strcmp(water_cond,'petzold_maalox')    
        
        load petzold_data_maalox_orig   
        %petzold_angle = [0.1:0.05:180]';
        angle_water = [0;angle_water];
        scatter_water = [scatter_water(1);scatter_water];
        cdf_scatter = cumtrapz(angle_water,scatter_water.*sind(angle_water));
        %cdf_scatter = interp1(angle_water,cdf_scatter,petzold_angle);
        cdf_scatter = cdf_scatter ./ max(cdf_scatter);
        petzold_angle_rad = angle_water.*pi./180;
        
    elseif strcmp(water_cond,'widemann_maalox')
        load widemann_maalox
        cdf_scatter = cumtrapz(widemann_maalox(:,1),widemann_maalox(:,2).*sin(widemann_maalox(:,1)));
        cdf_scatter = cdf_scatter ./ max(cdf_scatter);
        petzold_angle_rad = widemann_maalox(:,1);
    
        
    elseif strcmp(water_cond,'mie_1_micron')
        if (isunix()) 
            load('Berrocal/mieVSF1micron') 
        else
            load('Berrocal\mieVSF1micron') 
        end
        sphereAngles = sphereAngles.*pi./180;
        cdf_scatter = cumtrapz(sphereAngles,sin(sphereAngles).*mieVSF1micron);
        cdf_scatter = cdf_scatter ./ max(cdf_scatter);
        petzold_angle_rad = sphereAngles;
        
        
    elseif strcmp(water_cond,'petzold_harbor')
        load petzold_ocean
        cdf_scatter = cumtrapz(petzold_ocean(:,1),petzold_ocean(:,2).*sin(petzold_ocean(:,1)));
        cdf_scatter = cdf_scatter ./ max(cdf_scatter);
        petzold_angle_rad = petzold_ocean(:,1);
        
    elseif strcmp(water_cond,'petzold_coastal')
        load petzold_ocean
        cdf_scatter = cumtrapz(petzold_ocean(:,1),petzold_ocean(:,3).*sin(petzold_ocean(:,1)));
        cdf_scatter = cdf_scatter ./ max(cdf_scatter);
        petzold_angle_rad = petzold_ocean(:,1);
        
    elseif strcmp(water_cond,'petzold_clear')
        load petzold_ocean
        cdf_scatter = cumtrapz(petzold_ocean(:,1),petzold_ocean(:,4).*sin(petzold_ocean(:,1)));
        cdf_scatter = cdf_scatter ./ max(cdf_scatter);
        petzold_angle_rad = petzold_ocean(:,1);
        
    end

elseif strcmp(func_type,'calc')
        
    theta = [[0:0.01:10] [10.1:0.1:180]];
    
    if strcmp(water_cond,'hg')
        
        % HG phase function
        %g = 0.924;
        gsqr = g^2;       
        angle = theta.*pi./180;

        intensity = ((1/(4*pi)).*(1 - gsqr)) ./ (1 + gsqr - 2.*g.*cos(angle)).^(3/2);
        cdf_scatter = cumtrapz(angle,intensity.*sin(angle));
        cdf_scatter = cdf_scatter ./ max(cdf_scatter);
        petzold_angle_rad = angle;
        
    else

        if strcmp(water_cond,'harbor')
            % harbor water
            c = 2.190;
            a = 0.366;
        elseif strcmp(water_cond,'coastal')
            % coastal water
            c = 0.22 + 0.179;
            a = 0.179;
        else
            % Clear water
            c = 0.0374 + 0.114;
            a = 0.114;
        end

        b = c-a;
        albedo = (c-a)/c;   % Water albedo is scattering coef./atten. coef. (b/c unitless)


        q = 2.598 + 17.748*sqrt(b) - 16.722*b + 5.932*b*sqrt(b);
        k1 = 1.188 - 0.688*albedo;
        k2 = 0.1*(3.07 - 1.90*albedo);
        k3 = 0.01*(4.58 - 3.02*albedo);
        k4 = 0.001*(3.24 - 2.25*albedo);
        k5 = 0.0001*(0.84 - 0.61*albedo);

        phase_func = (4*pi/b)*exp(q*(1 + ((-1)^1*k1*theta.^(1/2)) ...
                                        + ((-1)^2*k2*theta.^(2/2)) ...
                                        + ((-1)^3*k3*theta.^(3/2)) ...
                                        + ((-1)^4*k4*theta.^(4/2)) ...
                                        + ((-1)^5*k5*theta.^(5/2))));
        %phase_func_adj = phase_func.*sind(theta);

        cdf_phase_func = cumtrapz(theta,phase_func.*sind(theta));
        cdf_phase_func = cdf_phase_func ./ max(cdf_phase_func);
        cdf_phase_func_interp = interp1(cdf_phase_func,theta,[0:0.001:1]);

        cdf_scatter = cdf_phase_func;
        petzold_angle_rad = theta.*pi./180;

    end
        
end
        
       

