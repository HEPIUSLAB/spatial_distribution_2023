% Function: ROI_from_manual
%
% Purpose: to fill in dorsal and ventral borders of the cord and create ROI
%
% Input parameters: 
%       vid: double (m x n x totalframes)
%       frames: double (1 x selectedframes)
%       dorsal: double (npoints x 2 x totalframes)
%       vid: double (npoints x 2 x totalframes)
%
% Output parameters:
%       spineROI: double/logical (m x n)
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)

function spineROI = ROI_from_manual(vid, dorsal, ventral, frames)

    if nargin == 3
        if length(size(vid))==3
            frames = 1:size(vid, 3);
        elseif length(size(vid))==4
            frames = 1:size(vid, 4);
        else
            error('Input is not a video');
        end
    end
    
    % initialize binary ROI image
    spineROI = zeros([size(vid, 1, 2), length(frames)]);

    for tt = 1:length(frames)        
        Xdors = [];
        Ydors = [];
        Xvent = [];
        Yvent = [];
        
        % convert dorsal Line ROIs into pixel indices
        for point = 1:size(dorsal, 1)-1
            [X, Y] = bresenham(dorsal(point,1,tt), dorsal(point,2,tt), dorsal(point+1,1,tt), dorsal(point+1,2,tt));

            Xdors = [Xdors; X];
            Ydors = [Ydors; Y];
        end
        
        % convert ventral Line ROIs into pixel indices
        for point = 1:size(ventral, 1)-1

            [X, Y] = bresenham(ventral(point,1,tt), ventral(point,2,tt), ventral(point+1,1,tt), ventral(point+1,2,tt));

            Xvent = [Xvent; X];
            Yvent = [Yvent; Y];
        end
        
        % fill in pixels between borders
        for xd = 1:length(Xdors)
            xv = find(Xvent == Xdors(xd),1,'first');

            if ~isempty(xv)
                spineROI(Ydors(xd):Yvent(xv), Xdors(xd), tt) = 1;
            end
        end

    end

end