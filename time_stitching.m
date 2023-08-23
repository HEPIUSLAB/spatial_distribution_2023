% Function: time_stitching
%
% Purpose: to stitch flow parameter data and link through time
%
% Input parameters: 
%       folder: string
%       surgery_data: structure (from loaded mat file)
%       paramname: string
%
% Output parameters:
%       t_USpre: double array (1xnumframes)
%       USpre: double array (1xnumframes)
%       t_USpost: double array (1xnumframes)
%       USpost: double array (1xnumframes)
%       inj_time
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)

function [t_USpre, USpre, t_USpost, USpost, inj_time] = time_stitching(folder, surgery_data, paramname)
    
    c = dir(folder);
    dsets = {c([c.bytes]~=0).name};
    
    % load old dicom listing
    listing = load(fullfile(folder,'..','..','DICOM listing.mat')).listing;
    listing(1).FlowTrace = [];
    listing(1).FlowTraceRed = [];

%     % error if parameter doesn't exist
%     if ~isfield(listing, paramname)
%         error('This parameter has not been calculated for this file')
%     end
    
    % load surgery date for record-keeping
    surg_date = listing(1).AcquisitionDate(3:end);
    if strcmp(surg_date,'221007')
        surg_date = '221006';
    end
    
    % for command window update
    disp('Loading ultrasound data...');
    lineLength = 0;
    
    %%%%% stitch data from same original DICOM back together %%%%%
    for ii = 1:length(dsets)
        
        % update command window
        if mod(ii, 25) == 0 || ii == length(dsets)
            fprintf(repmat('\b',1,lineLength));
            lineLength = fprintf('Loading File: %d of %d\n', ii, length(dsets));
        end
        
        % find original dataset for time information
        dataset = erase(erase(dsets{ii}, '.mat'), 'Velocity Map - ');    
        dicomind = find(strcmp({listing.name}, dataset(1:5)));
        
        
        % load flow parameter
        trace = load(fullfile(folder,sprintf('Velocity Map - %s.mat', dataset))).(paramname);
        
        % store flow parameter in listing
        listing(dicomind).FlowTrace = [listing(dicomind).FlowTrace; trace];

        close all    
    end
    
    
    %%%%%% REMOVING ERROR VIDEOS (from metadata) %%%%%%%
    vids_to_delete = surgery_data(strcmp({surgery_data.date}, surg_date)).vids_to_delete;    
    disp('Removing error DICOMs...');
    
    % erroneous first DICOM recorded on 221208
    if strcmp(surg_date,'221208')
        for ii = 1:length(listing)
            listing(ii).time = listing(ii).time + listing(1).NumberOfFrames*listing(1).FrameTime/1000;
        end
    end
    
    for ii = 1:length(vids_to_delete)
        listing(strcmp({listing.name}, vids_to_delete(ii))) = [];
    end
    
    
    %%%%%% STITCHING DICOMs and creating time axis %%%%%%%    
    disp('Stitching the DICOMs...')
    
    t_USpre = [];
    USpre = [];
    t_USpost = [];
    USpost = [];
    
    % determine injury DICOM index (from metadata)
    inj_ind = find(strcmp({listing.name}, surgery_data(strcmp({surgery_data.date}, surg_date)).inj_vid));
    
    % when no injury, processing is slightly different
    if isempty(inj_ind)
        inj_ind = length(listing);
        skip_var = 0;
    else
        skip_var = 1;
    end
    
    % calculate time of injury
    try
        inj_time = listing(inj_ind).time - listing(inj_ind).FrameTime*listing(inj_ind).NumberOfFrames/1000 ... 
                   + listing(1).FrameTime*listing(1).NumberOfFrames/1000;
    catch
        inj_time = 0;
    end
    
    
    % combine ultrasound for pre-injury
    for dicomind = 1:inj_ind-1   
        
        % SMI videos have ~20 frames of initialization (no signal)
        if strcmp(listing(dicomind).Mode,'SMI')
            cutoff = 20;
        else
            cutoff = 1;
        end
        
        % create time and flow variables
        t_USpre = [t_USpre ((cutoff-1):size(listing(dicomind).FlowTrace,1)-1) * listing(dicomind).FrameTime/1000 + ...
                listing(dicomind).time - listing(dicomind).NumberOfFrames*listing(dicomind).FrameTime/1000 + listing(1).FrameTime*listing(1).NumberOfFrames/1000];
        USpre = [USpre; listing(dicomind).FlowTrace(cutoff:end,:)];
    
    end
    
    % combine ultrasound for post injury
    for dicomind = inj_ind+skip_var:length(listing)   % for injury experiments
        
        % SMI videos have ~20 frames of initialization (no signal)
        if strcmp(listing(dicomind).Mode,'SMI')
            cutoff = 20;
        else
            cutoff = 1;
        end
        
        % create time and flow variables
        t_USpost = [t_USpost ((cutoff-1):size(listing(dicomind).FlowTrace,1)-1) * listing(dicomind).FrameTime/1000 + ...
                listing(dicomind).time - listing(dicomind).NumberOfFrames*listing(dicomind).FrameTime/1000 + listing(1).FrameTime*listing(1).NumberOfFrames/1000];
        USpost = [USpost; listing(dicomind).FlowTrace(cutoff:end,:)];
    
    end
end
