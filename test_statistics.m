% Old helper function for testing statistical calculations on the samples.

samples = 2.5*randn(1e6,1)+10;
count = 0;
meanval = 0;
sqravg =0;
variance = 0;

for i=1:length(samples)
   count =  count+1;
   meanval = ((count-1)*meanval + samples(i)) / count;   % E[theta] = (N*E[theta] + theta) / N+1 
   sqravg = ((count-1)*sqravg + samples(i)^2)/count;     % E[theta^2]
   variance = (count/(count-1))*(sqravg - meanval^2);    % Var[theta] = E[theta^2] - E[theta]^2
   
   
end




actual_var = var(samples)
variance
actual_mean = mean(samples)
meanval