% beam.m

% A MATLAB program
%   Steven L. Jacques, Dec. 1998.
%
% Generates the probability density function p(r)
% and the relative irradiance E(r) [mm^-2]
% for a 2D Gaussian beam profile associated with
% a laser beam. Histograms based on 10,000 random numbers
% are created, and compared with the analytic expressions.
%
% Two figures are generated: p(r) and E(r).
%

close
close
clear

b = 1;                       % arbitrary choice of 1/e beam radius. 

N = 10000;                   % number of random numbers used
rnd1 = rand(1,N);            % create array of random numbers
r1 = b*sqrt(-log(1 - rnd1)); % array of selected values r1
                             % Note: log() is base e

nb = 30;                     % number of bins
[n,x] = hist(r1,nb);         % n(x) = histo(x=r1) using nb bins
dx = x(3)-x(2);              % bin size
n = n/N/dx;                  % normalize by N and dx
sum(n*dx)                    % verify that integral equals unity
nr = n./(2*pi*x);            % normalize by circumference of ring at r=x

r = (0:.1*b:3*b);            % create an array of r values
E = exp(-r.^2/b^2)/(pi*b^2); % calculate an array of E values,
                             %   using a pure 2D Gaussian.

% Create figure of p(r)
figure
plot(x, n,'yo')
hold on
plot(r, E*2*pi.*r,'r-')
xlabel('r')
ylabel('p(r)')
xlabel('r [mm]')
title(['b = ' int2str(b) ' [mm]'])
hold off

% Create figure of E(r)
figure
plot(x,nr,'o')
xlabel('r')
ylabel('p(r)')
hold on
plot(r,E,'r-')
title(['b = ' int2str(b) ' [mm]'])
ylabel('E(r) [mm^-2]')
xlabel('r [mm]')
hold off
