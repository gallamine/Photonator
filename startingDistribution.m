% Starting intensity distribution for the data from E. Barrocal's thesis
% and publications. Used to verify the MC simulation output

function [xPos,yPos] = startingDistribution(numPhotons)

    minVal = 1;
    maxVal = 200;
    
    %im = imread('Berrocal/laser_intensity.bmp');
    if (isunix())
        load('Berrocal/Data/1micOD2/Incident');
    else
        load('Berrocal\Data\1micOD2\Incident');
    end
    Incident(Incident == 0) = 1;

    
%     Incident = double(rgb2ind(im, jet,'nodither'));
%     Incident(Incident==0) = 0.00001*rand();
 
    imSize = size(Incident);

    laserIntensityUnroll = reshape(Incident,imSize(1)*imSize(2),1);

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
       xPos(j) = ceil(vectPos(j)/imSize(1));
       yPos(j) = floor(mod(vectPos(j)-1,imSize(2))) + 1;
       if xPos(j) == 201
           disp('oh!')
       end
%        
%        if  ((xPos(j) > 200) || (yPos(j) > 200))
%            disp('Woah nelly!')
%        end 
    end

    xPos = (xPos - 100.5) ./ 99.5;      % scale to +-1
    yPos = (yPos - 100.5) ./ 99.5;      % scale to +-1

end