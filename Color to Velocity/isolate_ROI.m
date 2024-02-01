clear; close all;
% Purpose: Create subROIs (cranial, injury, and caudal) to run through
% Scholbach Scholbach analysis
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)


%% Section 1: Draw the ROIs

% root_path = uigetdir;
% root_path = 'G:\Data\Rat SCBF AH';
% root_path = 'F:\Data\220930 Rat SCAR w propofol';
% root_path =  'F:\Data\221117 Pig SCAR';
% root_path =  'E:\Data\221208 Rat SCAR with epi';
root_path =  'D:\Data\230322 Pig SCAR';


a = dir(root_path);
expsets = {a(~contains({a.name}, '.')).name};

for ii = 1:length(expsets)
    folder = fullfile(root_path, expsets{ii});

    % skips folders without listing
    try
        full_dicom_listing = load(fullfile(folder, 'DICOM listing.mat')).listing;
    catch
        continue;
    end
    
    folder = fullfile(folder, full_dicom_listing(1).AcquisitionDate);
    
%     if exist(fullfile(folder, 'multiROI.mat'), "file")
%         continue;
%     end

    
    drawROIs(full_dicom_listing, folder);
end


%% Section 2: Process through Scholbach Scholbach

% root_path = uigetdir;
% root_path = 'G:\sample';

a = dir(root_path);
expsets = {a(~contains({a.name}, '.')).name};

for ii = 1:length(expsets)
    folder = fullfile(root_path, expsets{ii});

    % skips folders without listing
    try
        full_dicom_listing = load(fullfile(folder, 'DICOM listing.mat')).listing;
    catch
        continue;
    end
    
    folder = fullfile(folder, full_dicom_listing(1).AcquisitionDate);
    
    run_perfusion(full_dicom_listing, folder);
    rmdir(fullfile(folder, 'Quantitative Data', 'Injury', full_dicom_listing(1).AcquisitionDate),'s');
    
end



%% DO NOT RUN: function definitions and testing
% 
% folder = 'C:\Users\Denis\Documents\JHSOM\PhD\Data\220126 Test Data for MATLAB\20211111';
% 
% full_dicom_listing = load(fullfile(folder,'..', 'DICOM listing.mat')).listing;
% 
% drawROIs(full_dicom_listing, folder);
% run_perfusion(full_dicom_listing, folder);
% rmdir(fullfile(folder, 'Quantitative Data', 'Injury', full_dicom_listing(1).AcquisitionDate),'s');
% rmdir(fullfile(folder, 'Quantitative Data', 'Cranial', full_dicom_listing(1).AcquisitionDate),'s');
% rmdir(fullfile(folder, 'Quantitative Data', 'Caudal', full_dicom_listing(1).AcquisitionDate),'s');

function drawROIs(full_dicom_listing,folder)
    
    vid_inds = find([full_dicom_listing.NumberOfFrames]>1);
    
    if ~exist(fullfile(folder, 'multiROI'), 'file')
        save(fullfile(folder, 'multiROI'),"vid_inds");
    else
        load(fullfile(folder, 'multiROI'));
    end

    if exist('ROI_inj', 'var')
        ROI_inj = struct();
    end

    for vid = 1:length(vid_inds)
        img = full_dicom_listing(vid_inds(vid)).image;
        
        if vid == 1 | vid == 22         % for one ROI for all vids
%         if true             % for one ROI for each DICOM
            if exist('multiROI.png', "file")        
                f = figure;
                imshow(imread('multiROI.png'));
            end
        
            f1 = figure('visible','on'); % Display image
            imagesc(img);axis off;
            ax = gca; % Identify current axes
        %     title(dataset);
            
            % Confirm user satisfaction with drawn ROI and redraw if not good
            user_satisfied = false;
            while ~user_satisfied
                hinj = drawfreehand(ax); % Draw ROI on the image    
                answer = questdlg('Are you satisfied with the injury ROI?');
            
                switch answer
                    case 'Yes'
                        user_satisfied = true;
                    case 'No'
                        delete(hinj);
                    case 'Cancel'
                        close(f1);
                        throw(MException('roiDraw:userCancel', 'User canceled run'));
                end
            end
            
            
            
            row_borders_inj = round(hinj.Position(:,2)); % Isolate the row borders of the ROI
            col_borders_inj = round(hinj.Position(:,1)); % Isolate the column borders of the ROI
            
            
            saveas(f1, 'multiROI.png');
            close(f1);
        
            try close(f)
            end
        end


        ROI_inj(vid).ROI = borders2inds(row_borders_inj, col_borders_inj);
           
        if vid == length(vid_inds)
            save(fullfile(folder, 'multiROI'), 'ROI_inj','-append');
        end
    end
end


function run_perfusion(full_dicom_listing, folder)
    load(fullfile(folder, 'multiROI.mat'));
    a = dir(fullfile(folder, 'Velocity Maps'));
    dsets = {a([a.bytes]~=0).name};
    date = full_dicom_listing(1).AcquisitionDate;
    
    for ii = 1:length(dsets)
        
        dataset = erase(erase(dsets{ii}, '.mat'), 'Velocity Map - ');
    %     disp(dataset);
        
        vid = find(strcmp({full_dicom_listing.name}, dataset(1:5)));
        
        if ~ismember(vid, vid_inds)
            continue
        else
            vid = find(vid_inds==vid);
        end

%         disp(vid)



        if exist(fullfile(folder, 'Quantitative Data','Injury', sprintf('Injury Quant - %s.mat',dataset)),'file')
            continue;
        end
        
        % Load calculated velocity map from 'main' script
        load(fullfile(folder,'Velocity Maps',sprintf('Velocity Map - %s.mat',dataset)));
        
        
        fprintf('Velocity Map - %s - Injury\n',dataset);
        scanner_parameters.folder = fullfile(folder,'Quantitative Data', 'Injury');
        if ~exist(fullfile(scanner_parameters.folder,date, 'Quantitative Data'),'dir')
            mkdir(fullfile(scanner_parameters.folder,date, 'Quantitative Data'));
        end
        perfusion_display_function(scanner_parameters, fliplr(ROI_inj(vid).ROI), velocity);
        movefile(fullfile(scanner_parameters.folder,date, 'Quantitative Data',sprintf('Matlab Quant - %s.mat',scanner_parameters.curr_dataset)), ...
            fullfile(scanner_parameters.folder,sprintf('Injury Quant - %s.mat',scanner_parameters.curr_dataset)));
        

    
    end
end


function ROI_idx = borders2inds(lateral_borders, axial_borders)
    
    new_borders = []; %Initialize matrix for ROI boundary
    for i = 1:length(lateral_borders)-1
        x1 = lateral_borders(i); % Get [x,y] coordinates for current row and the following row
        y1 = axial_borders(i);
        x2 = lateral_borders(i+1);
        y2 = axial_borders(i+1);
    
        [x,y] = bresenham(x1,y1,x2,y2); % Use bresenham algorithm to fill in pixel values between two points on the line in the ROI boundary
        new_borders = [new_borders; x y]; % Add new points to the new_borders matrix
    end
    % connect end point to starting point
    x1 = lateral_borders(end);
    y1 = axial_borders(end);
    x2 = lateral_borders(1);
    y2 = axial_borders(1);
    [x,y] = bresenham(x1,y1,x2,y2);
    new_borders = [new_borders; x y];
    
    lat_new = new_borders(:,1); % Separate lateral and axial borders
    ax_new = new_borders(:,2);
    
    prev_axial_range = [];  %Initialize matrix
    
    ROI_idx = []; %Initialize ROI indices
    for x = min(lat_new): max(lat_new) % Loop through all of the columns contained in the ROI
        idx = find(lat_new==x); % Find all ROI indices associated with current column
        if ~isempty(idx) % If there are ROI indices associated with the current column, continue
            axial_idx = ax_new(idx); % Pull the row indices associated with the current column index
            axial_max = max(axial_idx); % Find the last row in the ROI for the current column
            axial_min = min(axial_idx); % Find the first row in the ROI for the current column
            axial_range = axial_min:axial_max; % Define the range of rows in the ROI for the current column
            if size(axial_range,2)>=1 % If there is more than one pixel in this range of rows, continue
                ROI_idx = [ROI_idx; repmat(x,size(axial_range')) axial_range']; % Add current range of rows to the ROI_idx matrix
            end
        end
    end
end




