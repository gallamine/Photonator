% Empirical Petzold phase function from:
% THEORETICAL AND EMPIRICAL PHASE FUNCTIONS FOR MONTE CARLO CALCULATIONS OF
% LIGHT SCATTERING IN SEAWATER by Vladimir I. Haltrin

theta = [[0:0.01:10] [10.1:0.1:180]];
% theta = [0:180/5000:180];
% theta = [0:pi/5000:pi];
 

% % harbor water
% c = 2.190;
% a = 0.366;
% % coastal water
% c = 0.22 + 0.179;
% a = 0.179;
% % Clear water
% c = 0.0374 + 0.114;
% a = 0.114;
% Petzold Clear 1
c = 0.199;
a = 0.082;
% % Petzold Coastal 1
% c = 0.470;
% a = 0.195;


b = c-a;
albedo = (c-a)/c   % Water albedo is scattering coef./atten. coef. (b/c unitless)


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
                            
beta = (b/(4*pi)).*phase_func;
%phase_func_adj = phase_func.*sind(theta);

%%

normalization_eq = 0.5.*trapz(theta,(phase_func.*sin(theta*pi/180)))
% normalization_eq = 0.5.*trapz(theta,(phase_func.*sind(theta)))
  
%%

cdf_phase_func = cumtrapz(theta,phase_func.*sind(theta));
cdf_phase_func = cdf_phase_func ./ max(cdf_phase_func);
cdf_phase_func_interp = interp1(cdf_phase_func,theta,[0:0.001:1]);

% hold on;                            
% loglog(theta,phase_func)