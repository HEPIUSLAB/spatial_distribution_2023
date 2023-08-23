% Function: flow_pixelwise
%
% Purpose: compile the pixelwise analysis in the structures. Written after
% I had already made some really messy data structures
%
% Input parameters:
%       S_NC: non-contrast structure
%       S_FM: FlowMorph structure (has injury severity)
%       prepost: whether its a pre or post structure
%       plot_bool: logical
%
% Output parameters:
%       S_NC: updated non-contrast structure
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)

function S_NC = flow_pixelwise(S_NC, S_FM, prepost, plot_bool)

if nargin == 2
    plot_bool = false;
end

for ii = 1:length(S_NC)
        
        
    % define injury as line drawn using draw_injury.m
    if strcmp(prepost, 'pre')
        injLocRound = round(S_FM(ii).InjuryLocationPre-S_FM(ii).ROIWindowPre(1,:));
    else
        injLocRound = round(S_FM(ii).InjuryLocationPost-S_FM(ii).ROIWindowPost(1,:));
    end
    
    % load necessary variables
    delxy = S_FM(ii).delxy;
    NCUSim = S_NC(ii).velim;
    spineROI = S_NC(ii).SpineROI;
    NCUSim(~spineROI) = 0;

    % distance function per pixel
    injIm = zeros(size(NCUSim));
    
    % this try deals with failures due to ROI formation
    try
        injIm(sub2ind(size(injIm),injLocRound(:,1),injLocRound(:,2))) = 1;        
    catch
        injIm(sub2ind(size(injIm),injLocRound(:,1)+20,injLocRound(:,2)+20)) = 1;
    end
    
    % plot if its asked for
    if plot_bool
        subplot(ceil(length(S_FM)/3),3,ii)
        imshow(imdilate(injIm, strel('disk', 10))+NCUSim)
    end
    
    % calculate distance image
    distIm = bwdist(injIm)*delxy*10;
    distIm(:,1:floor(mean(injLocRound(:,2)))) = -distIm(:,1:floor(mean(injLocRound(:,2))));
    distIm(~spineROI) = nan;

    % ignore pixels not in ROI
    distNCUS = distIm(NCUSim>0);
    
    % measured flow of each pixel in ROI
    NCUSflo = NCUSim(NCUSim>0);
    
    % store variables back in structure
    S_NC(ii).InjuryImage = injIm;
    S_NC(ii).DistanceImage = distIm;
    S_NC(ii).DistanceNCUS = distNCUS;
    S_NC(ii).NCUSFlow = NCUSflo;
            
end
end