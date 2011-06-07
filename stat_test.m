% Helper program to determine the validity of various ways for calculating
% the mean and variance of the data. Online, offline, and grouped methods.

%%
num_photons = 3000;
num_rx = 1000;
numTerm = num_photons - num_rx;
values = [rand(num_rx,1);zeros(numTerm,1)];

mean1 = mean(values);
var1 = var(values);


% % Partition values into sets of equal size
% 
% var3 = var(values(1:num_rx));
% g = 10;
% k = num_rx / g;
% means = mean(reshape(values(1:num_rx),k,g));
% vars = var(reshape(values(1:num_rx),k,g));
% 
% var3est = ((k - 1)/(g*k-1))*(sum(vars) + k*(g-1)*var(means)/(k-1)); 
% 
% var3 - var3est

% % Partition values into sets of unequal size
% 
% g = 5;
% k = [200 200 300 300 2000];
% index = cumsum(k);
% 
% means(1) = mean(values(1:index(1)));
% vars(1) = var(values(1:index(1)));
% for j = 2:length(k)
%    means(j) = mean(values(index(j-1)+1:index(j)));
%    vars(j) = var(values(index(j-1)+1:index(j)));
% end
% 
% 
% var3est = (1/(num_photons-1)).*(sum((k-1).*vars) + sum(k.*(means - mean1).^2)); 
% 
% var1 - var3est

% Incremental calculation of variance/mean from partitions of unequal sizes

g = 10;
k = [300 300 300 300 300 300 300 300 300 300];

index = cumsum(k);

% mean/variance returned from function (partition mean and variances)
means(1) = mean(values(1:index(1)));        
vars(1) = var(values(1:index(1)));          
for j = 2:length(k)                             
   means(j) = mean(values(index(j-1)+1:index(j)));
   vars(j) = var(values(index(j-1)+1:index(j)));
end

%incremental calculation of mean and variance (arbitrarry partition)
totalMean = 0;
vSum = 0;
diffSum = 0;
for j = 1:length(k)
    delta = means(j) - totalMean;               % incremental calculation of total sample mean
    totalMean = totalMean + (k(j)*delta)/index(j);  % weight mean calculation
    vSum = vSum + (k(j)-1)*vars(j) + k(j)*(delta)*(means(j) - totalMean);
    %diffSum = diffSum;

end

var3est = (1/(num_photons-1)).*(vSum);
var1 - var3est

%incremental calculation of mean and variance (constant partition)
totalMean = 0;
vSum = 0;
diffSum = 0;
kCon = k(1);
for j = 1:length(k)
    delta = means(j) - totalMean;               % incremental calculation of total sample mean
    totalMean = totalMean + (delta)/j;  % weight mean calculation
    vSum = vSum + vars(j) + (kCon/(kCon-1))*(delta)*(means(j) - totalMean);
    %diffSum = diffSum;

end

var4est = ((kCon-1)/(num_photons-1)).*(vSum);
var1 - var4est

% % Partition of two, unequal, sizes (last partition, all zeros) - correct
% M2 = (num_rx-1)*var(values(1:num_rx))  + (mean(values(1:num_rx)))^2*((num_rx*numTerm)/(num_photons));
% var3 = M2/(num_photons-1);
% 
% var3 - var1;
% 
% % online variance - correct!
% 
% mean2 = 0;
% M2 = 0;
% 
% for i = 1:num_rx
%     delta = values(i) - mean2;
%     mean2 = mean2 + delta/i;
%     M2 = M2 + delta*(values(i) - mean2);
% end
% 
% varRX = M2/(num_rx-1);
% varRX - var(values(1:num_rx));