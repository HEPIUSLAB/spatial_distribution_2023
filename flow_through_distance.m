% Function: flow_through_distance
%
% Purpose: calculate the flow at each distance (should really be named
% distance function but I don't feel like finding all the places I used it
%
% Input parameters: 
%       injROI: position of Line object
%       US: US sample image
%       delxy: distance per pixel (in cm)
%
% Output parameters:
%       distIm: image where pixel value is distance from injury (in mm)
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)

function distIm = flow_through_distance(injROI, US, delxy)
    
    injIm = zeros(size(US));
    
    % convert line position to pixel indices
    [injLocY,injLocX] = bresenham(injROI(1),injROI(3),injROI(2),injROI(4));
    injLocRound = round([injLocX,injLocY]);
    
    % create binary image of injury
    try
        injIm(sub2ind(size(injIm),injLocRound(:,1),injLocRound(:,2))) = 1;
    catch
        deleteinds = [find(injLocRound(:,1)<1); find(injLocRound(:,2)<1)];
        injLocRound(deleteinds,:)=[];
        injIm(sub2ind(size(injIm),injLocRound(:,1),injLocRound(:,2))) = 1;
    end
    
    % distance funciton of image
    distIm = bwdist(injIm)*delxy*10;
    distIm(:,1:round(mean(injLocRound(:,2)))) = -distIm(:,1:round(mean(injLocRound(:,2))));

end
