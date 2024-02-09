clear; close all;
% Purpose: Run only this script to run through entire processing pipeline
% automatically
% MAKE SURE TO SET ALL AUTOMODE VARIABLES TO TRUE!!! 
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)
% 

% for debugging
% fpath = 'C:\Users\Denis\Documents\Data\220124 single subject\';
% full_dicom_listing = extract_dicom_info(fpath, true);


%%%%%%%% MAKE SURE TO SET ALL AUTOMODE VARIABLES TO TRUE!!! %%%%%%%%%%%
%%%%%%%%%    DISREGARD IF YOU HAVE NOT EDITED THE CODE     %%%%%%%%%%%%

fpath = uigetdir();

% extract dicom frames and save autorun file
full_dicom_listing = extract_dicom_info(fpath, true);
% full_dicom_listing = extract_dicom_info('G:\20220329\170554',true);
folder = full_dicom_listing(1).folder;
listing = full_dicom_listing;
save(fullfile(folder, 'DICOM listing.mat'), 'listing');
save('folder_auto.mat','folder');

% preprocessing (this occurs before image saving) and update autorun folder
preprocessing
folder = fullfile(folder, full_dicom_listing(1).AcquisitionDate);
save('folder_auto.mat','folder');

% export images
export_dicom_frames(full_dicom_listing, 1, 'all');
% export_frames_fast(full_dicom_listing)

% run main and update autorun folder
main
save('folder_auto.mat','folder');

% post process to get Scholbach Scholbach numbers
postprocessing