% Function: color_mapping_function_large
%
% Purpose: split large DICOMs into processable chunks to input into
% color_mapping_function
%
% Input parameters:
%   foldername: char
%   listing: struct (directory files)
%       listing.name: char
%       listing.folder: char
%       listing.date: char
%       listing.bytes: double
%       listing.isdir: logical
%       listing.datenum: double
%   ROI_idx: double (matrix containing row and column indices for pixels in
%   the ROI)
%   img_size: struct
%       img_size.m_img: double (number of rows in the image)
%       img_size.n_img: double (number of columns in the image)
%   scanner_parameters: struct
%       scanner_parameters.curr_colorbarbounds: double
%       scanner_parameters.curr_dataset: char
%       scanner_parameters.curr_img_mode: char
%       scanner_parameters.curr_scanner: char
%       scanner_parameters.date: char
%       scanner_parameters.cmap: double
%
% Output parameters
%   none
%
% Last edited: 2/15/2022 (Denis Routkevitch, droutke1@jhmi.edu)
%       Added detection of completed analyses without overwrite option.
%       Possible addition in future

function color_mapping_function_large(foldername,listing,ROI_idx,img_size,scanner_parameters)

% split images into batches of 150
for ii = 3:150:length(listing)
    % determine output naming
    submap = (ii-3)/150+1;
    
    % check for already complete analyses (no overwriting in this function)
    if exist(fullfile(foldername,'..','..','Velocity Maps', ...
        sprintf('Velocity Map - %s - %04d.mat',scanner_parameters.curr_dataset, submap)), 'file')
        continue;
    end
    
    fprintf('Velocity Map - %s - %04d.mat\n',scanner_parameters.curr_dataset, submap);
    fprintf('Processing map %d/%d\n',submap, floor(length(listing)/150)+1);
    
    % determine images to analyze
    sublisting = listing([1,2,ii:min(length(listing),ii+149)]);
    imagenums = (ii:min(length(listing),ii+149))-2;

    tic
    if scanner_parameters.curr_img_mode == 'SMI' || scanner_parameters.curr_img_mode == 'CHI'
        [velocity,number_frames] = color_mapping_function_fast(foldername,sublisting,...
            ROI_idx,img_size,scanner_parameters); % Calculate velocity map for the current dataset
    else
        [velocity,number_frames] = color_mapping_function(foldername,sublisting,...
            ROI_idx,img_size,scanner_parameters); % Calculate velocity map for the current dataset
    end

    save(fullfile(foldername,'..','..','Velocity Maps',sprintf('Velocity Map - %s - %04d',scanner_parameters.curr_dataset, submap)), ...
        'velocity','number_frames','scanner_parameters','ROI_idx','-v7.3','imagenums'); % Save velocity map for the current dataset
    toc
    
end
end