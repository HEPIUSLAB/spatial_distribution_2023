% Function: draw_borders
%
% Purpose: draw borders through all different identified positions,
% starting with last. Use PolyLine
%
% Input parameters:
%   vid: double (mxnx3xnumframes or mxnxnumframes)
%   frames: double vector (the different positions)
%
% Output parameters:
%   Borders: double (set of indices)
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)

function Borders = draw_borders(vid, frames)

    Borders = [];

    if nargin == 1
        if length(size(vid))==3
            frames = 1:size(vid, 3);
        elseif length(size(vid))==4
            frames = 1:size(vid, 4);
        elseif length(size(vid))==2
            frames = 1;
        else
            error('Input is not an image or video');
        end
    end
    
    for ff = fliplr(frames)
        
        if length(size(vid))==3 || length(size(vid))==2
            I = vid(:,:,ff);
        elseif length(size(vid))==4
            I = vid(:,:,:,ff);
        else
            error('Input is not a video');
        end        
        
        % In last position, establish border
        if ff == frames(end)
            
            % last image
            f1 = figure;
            imshow(I);        
            ax = gca; % Identify current axes
            
            % draw border ROI
            user_satisfied = false;
            while ~user_satisfied
                h = drawpolyline(ax); % Draw ROI on the image    
                answer = questdlg('Are you satisfied with this ROI?');
            
                switch answer
                    case 'Yes'
                        user_satisfied = true;
                    case 'No'
                        delete(h);
                    case 'Cancel'
                        close(f1);
                        throw(MException('roiDraw:userCancel', 'User canceled run'));
                end
            end
            
            injpos= h.Position;
        
        % in all other positions, move border
        else        
            % copy previously drawn ROI onto pre image
            set(gcf, "Position", [1400 100 800 800])
            
            % show image
            f2 = figure;    
            imshow(I)
                
            ax = gca; % Identify current axes
            
            % Draw border ROI
            h = images.roi.Polyline(ax, 'Position', injpos); % Draw ROI on the image
            c = uicontrol('String','Continue','Callback',@button_pushed);
            uiwait(f2)    
            injpos= h.Position;
    
            close(f1)
            f1 = figure(f2);            
        end
        
        % store borders
        Borders = cat(3, injpos, Borders);
    end
    
    close(f1)
    
    % to exit drawing
    function button_pushed(src, event)
        uiresume(f2)
    end
end


