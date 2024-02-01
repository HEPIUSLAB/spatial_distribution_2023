clear; close all;
% Purpose: Complete post-processing (re-order velocity map, calculate mean
% velocity over the cardiac cycle, and create synchronized videos)
%
% Created by: Kelley Kempski (kkempski@jhmi.edu)
% Last edited: 1/25/22 (Denis Routkevitch droutke1@jhmi.edu)
%       Added looping and different folder selection mode
% Edit: 2/7/2022 added compatibility with DICOM splitting and new image
%       naming system
% Last edited: 2/15/2022 (Denis Routkevitch, droutke1@jhmi.edu)
%       Added folder change detection. Implemented new call to
%       perfusion_display_function (velocity and ROI_idx)
% Edit 4/5/2022  (DR):
%       Added autorun capability


%%%%
% for debugging convenience, create your own 'Sample DICOM Listing.mat' 
% by copying and renaming from export_dicom_frames() save
% during preprocessing
debugmode = false;


%%% SET THIS TO TRUE FOR AUTORUN %%%
automode = true;

% for overwrite
group = 'Save_Options';
try rmpref(group);
end

if debugmode
    % skips folder browser step
    full_dicom_listing = load('Sample DICOM Listing.mat').listing;
    folder = fullfile(full_dicom_listing(1).folder, full_dicom_listing(1).AcquisitionDate, 'Velocity Maps');
elseif automode
    % pulls from auto folder save
    load('folder_auto.mat', 'folder');
    folder = fullfile(folder, '..', 'Velocity Maps');
else
    [~, folder] = uigetfile('*.*');
end

% extract mat filenames
a = dir(folder);
dsets = {a([a.bytes]~=0).name};


% % Older version: UPDATE LINE 9: Surgery date (use same formatting as the folder name)
% date = '11-11-2021'; 

% % Older version: UPDATE LINE 12: Datasets you want to loop through
% for curr_dataset = 35
for ii = 1:length(dsets)
    
    dataset = erase(erase(dsets{ii}, '.mat'), 'Velocity Map - ');
%     disp(dataset);
    
    % check for existing mats and ask about overwrite option
    if exist(fullfile(folder,'..','Quantitative Data', ...
        sprintf('Matlab Quant - %s.mat',dataset)), 'file')
        
        if automode
            % automode does not overwrite
            continue;
        end
        
        pref =  'Overwrite';
        title = 'Overwrite?';
        
        quest = {sprintf('Post-processing for %s has already been done (although it may be incomplete).', ...
            dataset), 'Would you like to overwrite?'};
        pbtns = {'Yes','No'};
        
        pval = uigetpref(group,pref,title,quest,pbtns, ...
            'CheckboxString', 'Do this for all future datasets');
        
        switch pval
            case 'no'
                continue;
        end
    end
    
    
    % Load calculated velocity map from 'main' script
    load(fullfile(folder,sprintf('Velocity Map - %s.mat',dataset)));

    % skip calculation for CHI
    if scanner_parameters.curr_img_mode == 'CHI'
        continue;
    end
    
    % check for folder change and change necessary file paths if needed
    if ~strcmp(fullfile(folder,'..'), fullfile(scanner_parameters.folder,scanner_parameters.date, 'Velocity Maps','..'))
%         disp('hi')
        scanner_parameters.folder = fullfile(folder,'..','..');
    end
    
    % search for ROI_idx if not in loaded mat file
    if ~exist('ROI_idx', 'var')
        fullfile(folder, '..', 'PreProcessing',sprintf('PreProcessing - %s.mat', dataset(1:5)))
        load(fullfile(folder, '..', 'PreProcessing',sprintf('PreProcessing - %s.mat', dataset(1:5))), 'ROI_idx');
    end
    
    % Update scanner_parameters dataset name
    scanner_parameters.curr_dataset = dataset;
    % Save the new velocity map with updated scanner_parameters
    save(fullfile(folder,sprintf('Velocity Map - %s',dataset)), ...
        'velocity','number_frames','scanner_parameters','ROI_idx','-v7.3');
    

    % Save mean velocity and mean intensity for each frame
    perfusion_display_function(scanner_parameters, ROI_idx, velocity);


    % Create and save synchronized video for screenshots, calculated velocity
    % map, and velocity v. time plot
    % Takes 5-ever, comment out if you need
%     movie_saving_function(scanner_parameters)

    close all
end


% clean up workspace
try rmpref(group);
end
clear group pref title quest pbtns pval