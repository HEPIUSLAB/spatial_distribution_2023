clear; close all;
% Purpose: Run this script even more first (0/3) to extract pngs and
% necessary info from DICOM. Also, organize into folders for input into
% preprocessing script.
%
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)
% Last edited: 1/24/2022
%       Edited 1/28: commented out correct line (18 instead of 19)

% fpath = 'C:\Users\Denis\Documents\JHSOM\Long-Term Storage PhD\Data\220111 SCAR initial\2022_1_11_Ultrasound';
% fpath = 'C:\Users\Denis\OneDrive for Business\HEPIUS\Pig\Experiments\Terminal_Ultrasound_2021_12_16';
% fpath = 'C:\Users\Denis\Documents\JHSOM\Long-Term Storage PhD\Data\Test Data 2';
% fpath = 'C:\Users\Denis\OneDrive for Business\HEPIUS\Human OR\210126 Case 1\125741';
fpath = 'C:\Users\Denis\Documents\JHSOM\Long-Term Storage PhD\Data\220203 Pig SCAR';

% full_dicom_listing = extract_dicom_info(fpath, true);
full_dicom_listing = extract_dicom_info(true);

% export_dicom_frames(full_dicom_listing, 10, 1:3);
% export_dicom_frames(full_dicom_listing, 1, 'all');
export_frames_fast(full_dicom_listing)
