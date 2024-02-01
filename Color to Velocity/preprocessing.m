clear; close all;
% Purpose: Run this script first (1/3) to define scanner parameters and dataset
% information.
%
% Created by: Kelley Kempski (kkempski@jhmi.edu)
% Last edited: 1/25/2022 (Denis Routkevitch, droutke1@jhmi.edu)
%   Added compatibility with multi-DICOM processing
% Last edited: 2/15/2022 (Denis Routkevitch, droutke1@jhmi.edu)
%       Added folder change detection. Also added ability to recycle ROI
%       throughout entire recording session. (Autoreg experiments have
%       stable probe holder setup)
% Edit 4/5/2022  (DR):
%       Added autorun capability

%%%%
% for debugging convenience, create your own 'Sample DICOM Listing.mat' 
% by copying and renaming from export_dicom_frames() save
% during preprocessing
debugmode = false;  
readgui = false; % run code the original way (12/6/21)

oneROI = false; % use same ROI for each image
% firstDICOM = 1; % first dataset to start at, in case ROI changes during experiment

%%% SET THIS TO TRUE FOR AUTORUN %%%
automode = true;

% for overwrite
group = 'Save_Options';
try rmpref(group);
end

if readgui 
    % original parameter loading method (12/6/21)
    A = scanner_parameter_selection; % Trigger GUI for entering scanner and dataset information
    
    scanner_parameters = gui_read_function(A); % Read values from GUI

    scanner_parameters.folder = '..'; % back compatibility
    full_dicom_listing = 1;  % inserted for for loop functionality

elseif automode
    % pulls from auto folder save
    load('folder_auto.mat', 'folder');
    fpath = folder;
    full_dicom_listing = load(fullfile(folder, 'DICOM listing.mat')).listing;

elseif debugmode
    % skips folder browser step
    full_dicom_listing = load('Sample DICOM Listing.mat').listing;
else
    [fname, fpath] = uigetfile;
    full_dicom_listing = load(fullfile(fpath, fname)).listing;
end

firstDICOM = find((~strcmp([full_dicom_listing.Mode], 'BMODE') & [full_dicom_listing.NumberOfFrames]>1),1); % moved this to here for bug fix

for ii = firstDICOM:length(full_dicom_listing)
    
    % set bounds correctly for CHI (CEUS)
    if full_dicom_listing(ii).Mode == 'CHI'
        full_dicom_listing(ii).bounds = 1;
    end
    
    % since listing2params throws error if wrong modality, using try
    try
        if ~readgui
            % transcribe listing info to scanner_parameters (multi_DICOM
            % compatibility)
            scanner_parameters = listing2params(full_dicom_listing(ii));
        end
    catch
        % If not acceptable modality, preprocessing mat not created,
        % analysis not continued in further scripts
        continue;
    end
    
    % check for folder change and change necessary file paths if needed
    if ~strcmp(fullfile(fpath,'..'), fullfile(scanner_parameters.folder,'..'))
        scanner_parameters.folder = fpath;
    end

    % check for existing mats and ask about overwrite option
    if exist(fullfile(scanner_parameters.folder,scanner_parameters.date,'PreProcessing', ...
        sprintf('PreProcessing - %s.mat',scanner_parameters.curr_dataset)), 'file')
        
        
        pref =  'Overwrite';
        title = 'Overwrite?';
        quest = {sprintf('Pre-processing for %s has already been done.', ...
            scanner_parameters.curr_dataset), 'Would you like to overwrite?'};
        pbtns = {'Yes','No'};
        
        pval = uigetpref(group,pref,title,quest,pbtns, ...
            'CheckboxString', 'Do this for all future datasets');
        
        switch pval
            case 'no'
                load(fullfile(scanner_parameters.folder,scanner_parameters.date,'PreProcessing', ...
                    sprintf('PreProcessing - %s.mat',scanner_parameters.curr_dataset)), 'ROI_idx')

                continue;
        end
    end
    
    if automode
        % adds image to scanner parameters for setup_ROI and recycle_ROI
        scanner_parameters.img = full_dicom_listing(ii).image;
    end

    if oneROI && ii~=firstDICOM         % recycle ROI, extract other params
        [listing,img_size,foldername] = recycle_ROI(scanner_parameters);
    else                        % Draw ROI on single frame from the cine loop
        [ROI_idx,listing,img_size,foldername] = setup_ROI(scanner_parameters); 
    end
                                                                 
    save(fullfile(scanner_parameters.folder,scanner_parameters.date,'PreProcessing',sprintf('PreProcessing - %s',scanner_parameters.curr_dataset)), ...
        'scanner_parameters','listing','ROI_idx','img_size','foldername');

end

% clean up workspace
try rmpref(group);
end
clear group pref title quest pbtns pval
