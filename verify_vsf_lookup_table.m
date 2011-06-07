% Test the method for calculating scattering lengths and choosing values
% from the VSF (volume scattering function) distribution

num_photons = 1e6;
inv_b = 1;

[cdf_scatter,angle] = generate_scatter('measured','petzold_avg');

rand_array = rand(num_photons,3);    % generate a matrix for each photon with rand propogation, roll, and pitch
rand_array(:,3) = rand_array(:,3).*2.*pi;   %% Uniformly distributed over 0 -2Pi                                                                             




r = zeros(1,num_photons);
theta = zeros(1,num_photons);

% iterate over every single photon to calculate new position and
% whether it was received or not.
for i = 1:num_photons
        
    r(i) = -inv_b*log(rand_array(i,1));     % randomly generate optical path length      

    %Generate scattering angle from beam spread function
    k = 1;
    while rand_array(i,2) > cdf_scatter(k)
        k = k+1;
    end
    theta(i) = angle(k);
        
end

hold on;

figure(3);[N,X] = hist(r,100);
plot(X,N./N(1));
figure(4);[N,X] = hist(theta,500);
plot(X,N./N(1));
figure(5);[N,X] = hist(rand_array(:,3),360);
plot(X,N./N(1));