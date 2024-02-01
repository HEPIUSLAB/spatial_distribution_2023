% Function: setup_ROI
%
% Purpose: Open Dicom file for cine loops and save individual image frames. 
% Then, draw ROI for the current dataset and export values
%
% Input parameters:
%   scanner_parameters
%       scanner_parameters.curr_colorbarbounds: double
%       scanner_parameters.curr_dataset: char
%       scanner_parameters.curr_img_mode: char
%       scanner_parameters.curr_scanner: char
%       scanner_parameters.date: char
%       scanner_parameters.cmap: double
%
% Output parameters:
%   ROI_idx: double (matrix containing row and column indices for pixels in
%   the ROI)
%   listing: struct (directory files)
%       listing.name: char
%       listing.folder: char
%       listing.date: char
%       listing.bytes: double
%       listing.isdir: logical
%       listing.datenum: double
%   img_size: struct
%       img_size.m_img: double (number of rows in the image)
%       img_size.n_img: double (number of columns in the image)
%   foldername: char
%
% Created by: Kelley Kempski (kkempski@jhmi.edu)

% Edited: 3/23/2022 (Kelley Kempski)
%   Updated the drawn ROI boundaries by adding call to the bresenham
%   algorithm and updated ROI assignments accordingly

% Lasted edited: 1/25/22 (Denis Routkevitch, droutke1@jhmi.edu)
%       Added compatibility with any folder
%       Also added multi-processing functionality and confirm button
% Last edited: 2/15/2022 (Denis Routkevitch, droutke1@jhmi.edu)
%       Changed ROI drawing frame to last frame (dark frames issue)
% Edit 4/5/2022 (DR):
%       Added autorun capability


function [ROI_idx,listing,img_size,foldername] = setup_ROI(scanner_parameters)
date = scanner_parameters.date; % Extract surgery date
dataset = scanner_parameters.curr_dataset; % Extract dataset name
cmap = scanner_parameters.cmap; % Extract colormap
cmap_center = squeeze(cmap(:,9,:)); % Isolate center column of colormap for displaying velocity maps

% added option for selecting different folder
if isfield(scanner_parameters, 'folder')
    folder = scanner_parameters.folder;
else 
    folder = '..';
end

foldername     = fullfile(folder,date,'Images',dataset); % Set folder name for current dataset
listing        = dir(foldername); % Create directory of files contained in the folder
number_frames  = length(listing)-2; % Extract number of frames

% Load first frame of cine loop to select ROI.
if isfield(scanner_parameters, 'img')
    % for autorun, scanner parameters contains a frame
    img = scanner_parameters.img; 
else
    % for non autorun, scanner parameters does not contain a frame
    filename = fullfile(foldername,listing(end).name); 
    img = imread(filename); 
end

[m_img,n_img,~] = size(img); % Calculate size of the image
img_size.m_img = m_img; % Save number of rows to the img_size struct
img_size.n_img = n_img; % Save number of columns to the img_size struct

f1 = figure('visible','on'); % Display image
imagesc(img);axis off;
ax = gca; % Identify current axes
title(dataset);

% Confirm user satisfaction with drawn ROI and redraw if not good
user_satisfied = false;
while ~user_satisfied
    h = drawfreehand(ax); % Draw ROI on the image    
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

axial_borders = round(h.Position(:,2)); % Isolate the row borders of the ROI
lateral_borders = round(h.Position(:,1)); % Isolate the column borders of the ROI

close(f1);

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
if ~exist(fullfile(folder,date,'PreProcessing'),'dir')
    mkdir(fullfile(folder,date,'PreProcessing'));
end
if ~exist(fullfile(folder,date,'Velocity Maps'),'dir')
    mkdir(fullfile(folder,date,'Velocity Maps'));
end
if ~exist(fullfile(folder,date,'Quantitative Data'),'dir')
    mkdir(fullfile(folder,date,'Quantitative Data'));
end
if ~exist(fullfile(folder,date,'Videos'),'dir')
    mkdir(fullfile(folder,date,'Videos'))
end
end
