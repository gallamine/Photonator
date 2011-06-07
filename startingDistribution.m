% Starting intensity distribution for the data from E. Barrocal's thesis
% and publications. Used to verify the MC simulation output

function [xPos,yPos] = startingDistribution(numPhotons)

    minVal = 1;
    maxVal = 200;
    
    im = imread('Berrocal/laser_intensity.bmp');
    laserIntensity = double(rgb2ind(im, jet,'nodither'));
    laserIntensity(laserIntensity==0) = 0.00001*rand();

    imSize = size(laserIntensity);

    laserIntensityUnroll = reshape(laserIntensity,imSize(1)*imSize(2),1);

    imInt = cumtrapz(laserIntensityUnroll);

    imCDF = imInt./max(max(imInt));

    points = numPhotons;
    xPos = zeros(points,1);
    yPos = zeros(points,1);
    vectPos = zeros(points,1);
    for j = 1:points
        value = rand();
        %value = 0.838820707;
        closetIndex = lookupCDF(imCDF,value);
        vectPos(j) = (closetIndex-1) + (value - imCDF((closetIndex-1)))*(closetIndex - (closetIndex-1))/(imCDF(closetIndex) - imCDF(closetIndex-1));
       yPos(j) = ceil(vectPos(j)/imSize(1));
       xPos(j) = round(mod(vectPos(j)-1,imSize(2))) + 1;
%        
%        if  ((xPos(j) > 200) || (yPos(j) > 200))
%            disp('Woah nelly!')
%        end 
    end

    xPos = (xPos - 100.5) ./ 99.5;      % scale to +-1
    yPos = (yPos - 100.5) ./ 99.5;      % scale to +-1

end