% Function: recycle_ROI
%
% Purpose: When recycling ROI, extract all other necessary parameters in
% preprocessing
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
% Created by: Denis Routkevitch (droutke1@jhmi.edu)
%   Last edited: 2/9/2022
%       Basically scrounged the first part of setup_ROI for ROI extension
%       capability
% Edit 4/5/2022 (DR):
%       Added autorun capability


function [listing,img_size,foldername] = recycle_ROI(scanner_parameters)
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

% 
% % Load first frame of cine loop to select ROI
% filename = fullfile(foldername,listing(end).name);
% img = imread(filename); 

[m_img,n_img,~] = size(img); % Calculate size of the image
img_size.m_img = m_img; % Save number of rows to the img_size struct
img_size.n_img = n_img; % Save number of columns to the img_size struct
end