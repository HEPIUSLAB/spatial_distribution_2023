clear; close all;
% Purpose: Run only this script to run through entire processing pipeline
% automatically on multiple folders
% MAKE SURE TO SET ALL AUTOMODE VARIABLES TO TRUE!!! 
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)
% 

% for debugging
% fpath = 'C:\Users\Denis\Documents\Data\220124 single subject\';
% full_dicom_listing = extract_dicom_info(fpath, true);


%%%%%%%% MAKE SURE TO SET ALL AUTOMODE VARIABLES TO TRUE!!! %%%%%%%%%%%
%%%%%%%%%    DISREGARD IF YOU HAVE NOT EDITED THE CODE     %%%%%%%%%%%%

%% Section 1

root_path = uigetdir('D:\Data');
% root_path = 'G:\sample';

a = dir(root_path);
% expsets = {a(~contains({a.name}, '.')};
% expsets = {a(~contains({a.name}, '.')&contains({a.name}, 'US')).name};
expsets = {a(contains({a.name}, 'Elect')).name};
save('multi_auto.mat','expsets', "root_path");

for ii = 1:length(expsets)

    fpath = fullfile(root_path, expsets{ii});

    if isfile(fullfile(fpath, 'DICOM listing.mat'))
        continue;
    end

    % extract dicom frames and save autorun file
    full_dicom_listing = extract_dicom_info(fpath, true);
    folder = full_dicom_listing(1).folder;
    listing = full_dicom_listing;
    save(fullfile(fpath, 'DICOM listing.mat'), 'listing');

end

%% Section 2
% root_path = uigetdir('D:\Data');
% root_path = 'G:\sample';


% a = dir(root_path);
% % expsets = {a(~contains({a.name}, '.')};
% expsets = {a(~contains({a.name}, '.')&contains({a.name}, 'US')).name};
% save('multi_auto.mat','expsets', "root_path");

for ii = 1:length(expsets)
    save('ii_auto.mat','ii');
    folder = fullfile(root_path, expsets{ii});
    save('folder_auto.mat','folder');
    
    % preprocessing (this occurs before image saving) and update autorun folder
    preprocessing
    load('multi_auto.mat','expsets', "root_path");
    load('ii_auto.mat','ii');
    
end

%% Section 3

% root_path = uigetdir('D:\Data');
% root_path = 'G:\sample';
% root_path = 'E:\Data\230119 Rat SCBF and CEUS';

% a = dir(root_path);
% % expsets = {a(~contains({a.name}, '.')};
% expsets = {a(~contains({a.name}, '.')&contains({a.name}, 'US')).name};
% save('multi_auto.mat','expsets', "root_path");

for ii = 1:length(expsets)
    load('multi_auto.mat','expsets', "root_path");
%     folder = fullfile(root_path, expsets{length(expsets)-ii+1});
    folder = fullfile(root_path, expsets{ii});
    full_dicom_listing = load(fullfile(folder, 'DICOM listing.mat')).listing;
    folder = fullfile(folder, full_dicom_listing(1).AcquisitionDate);
    save('folder_auto.mat','folder');
    
    % export images
    export_dicom_frames(full_dicom_listing, 1, 'all');
    
    % run main and update autorun folder
    main
    save('folder_auto.mat','folder');
    
    % post process to get Scholbach Scholbach numbers
    postprocessing

    load('multi_auto.mat','expsets', "root_path");
    load('ii_auto.mat','ii');
end


%% Temporary garbage to save all snapshots in one place

% root_path = uigetdir;
root_path = 'G:\Data\Rat SCBF AH';

a = dir(root_path);
expsets = {a(~contains({a.name}, '.')).name};
% save('multi_auto.mat','expsets', "root_path");

for ii = 1:length(expsets)
%     load('multi_auto.mat','expsets', "root_path");
    disp(ii/length(expsets));
    folder = fullfile(root_path, expsets{ii});
    try
        full_dicom_listing = load(fullfile(folder, 'DICOM listing.mat')).listing;
    catch
        full_dicom_listing = extract_dicom_info(folder, true);
        listing = full_dicom_listing;
        save(fullfile(folder, 'DICOM listing.mat'), 'listing');        
    end
    
    inds = find([full_dicom_listing.Mode] == 'CDI'&[full_dicom_listing.NumberOfFrames]==1);

    if isempty(inds) || inds(1) == 1 || inds(1) == length(full_dicom_listing)
        inds = find([full_dicom_listing.Mode] == 'CDI');
        I = dicomread_debugged(fullfile(folder, full_dicom_listing(inds(end)).name),"Frames",full_dicom_listing(inds(end)).NumberOfFrames);
    else
        I = dicomread(fullfile(folder, full_dicom_listing(inds(1)).name));
    end

    imwrite(I, fullfile(root_path, 'All images', sprintf('%s.png', expsets{ii})));

%     folder = fullfile(folder, full_dicom_listing(1).AcquisitionDate);
%     save('folder_auto.mat','folder');
%     
%     % export images
%     export_dicom_frames(full_dicom_listing, 1, 'all');
%     
%     % run main and update autorun folder
%     main
%     save('folder_auto.mat','folder');
%     
%     % post process to get Scholbach Scholbach numbers
%     postprocessing
% 
%     load('multi_auto.mat','expsets', "root_path");
%     load('ii_auto.mat','ii');
end
