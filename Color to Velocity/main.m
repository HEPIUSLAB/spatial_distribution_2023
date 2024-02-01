clear; close all;
% Purpose: Create and save velocity maps
%
% Created by: Kelley Kempski (kkempski@jhmi.edu)
% Last edited: 1/25/22 (Denis Routkevitch droutke1@jhmi.edu)
%       Added looping and different folder selection mode
% Edit: 2/7/2022 added compatibility with DICOM splitting and new image
%       naming system
% Last edited: 2/15/2022 (Denis Routkevitch, droutke1@jhmi.edu)
%       Added folder change detection
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
%     full_dicom_listing = load('Sample DICOM Listing.mat').listing;
%     folder = fullfile(full_dicom_listing(1).folder, full_dicom_listing(1).AcquisitionDate, 'PreProcessing');

    folder = 'C:\Users\Denis\Documents\Data\220124 single subject\20210930\PreProcessing';
elseif automode
    % pulls from auto folder save
    load('folder_auto.mat', 'folder');
    folder = fullfile(folder, 'PreProcessing');
else
    [~, folder] = uigetfile('*.*');
end

% extract mat filenames
a = dir(folder);
dsets = {a([a.bytes]~=0).name};

% % Older version: UPDATE LINE 7: Surgery date (use same formatting as the folder name)
% surgery_date = '11-11-2021'; 

% % Older version: UPDATE LINE 10: Datasets you want to loop through

for ii = 1:length(dsets)
    
    dataset = dsets{ii};
    load(fullfile(folder, dataset)); % Load pre-processing information
    
    % check for folder change and change necessary file paths if needed
    if ~strcmp(fullfile(folder,'..'), fullfile(scanner_parameters.folder,scanner_parameters.date, 'PreProcessing','..'))
        scanner_parameters.folder = fullfile(folder,'..','..');
        foldername = fullfile(folder,'..','Images',scanner_parameters.curr_dataset);
    else
        foldername = fullfile(folder,'..','Images',scanner_parameters.curr_dataset);
    end
    
    
    % check for existing mats and ask about overwrite option
    if exist(fullfile(folder,'..','Velocity Maps', ...
        sprintf('Velocity Map - %s.mat',scanner_parameters.curr_dataset)), 'file')

        if automode
            % automode does not overwrite
            continue;
        end
        
        pref =  'Overwrite';
        title = 'Overwrite?';
        
        quest = {sprintf('Velocity mapping for %s has already been done.', ...
            scanner_parameters.curr_dataset), 'Would you like to overwrite?'};
        pbtns = {'Yes','No'};
        
        pval = uigetpref(group,pref,title,quest,pbtns, ...
            'CheckboxString', 'Do this for all future datasets');
        
        switch pval
            case 'no'
                continue;
        end
    end


    % no DICOM splitting (uncomment this and comment below section)
    tic
    % be mindful of fast vs original function
%     [velocity,number_frames] = color_mapping_function(foldername,dir(foldername),...
%         ROI_idx,img_size,scanner_parameters); % Calculate velocity map for the current dataset
%     
% 
%     save(fullfile(folder,'..','Velocity Maps',sprintf('Velocity Map - %s',scanner_parameters.curr_dataset)), ...
%         'velocity','number_frames','scanner_parameters','ROI_idx','-v7.3'); % Save velocity map for the current dataset
    

    % DICOM splitting (uncomment this and comment above section)
    color_mapping_function_large(foldername,dir(foldername),ROI_idx,img_size,scanner_parameters);

    toc
end


% clean up workspace
try rmpref(group);
end
clear group pref title quest pbtns pval