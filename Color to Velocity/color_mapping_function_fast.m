% Function: color_mapping_function_fast
%
% Purpose: Create 3D matrix of the velocity map for the current dataset.
%       Uses faster linear regression model
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
%   velocity: double (matrix containing 3D map of velocities - m x n x r
%   where m and n are the rows and columns, respectively, and r is the
%   number of frames
%   number_frames: double
%
% Created: 4/55/2022 (Denis Routkevitch, droutke1@jhmi.edu)
%       vectorization and added SMI linear model.

%%%%%% CURRENTLY ONLY COMPATIBLE WITH SMI!!!!! %%%%%%%

function [velocity,number_frames] = color_mapping_function_fast(foldername,listing,ROI_idx,img_size,scanner_parameters)

% temporary warning if not SMI
if ~strcmp(scanner_parameters.curr_img_mode, 'SMI') && ~strcmp(scanner_parameters.curr_img_mode, 'CHI')
    warning(['This set (%s) is not SMI, which is not yet compatible' ...
        ' with the fast function. Velocity map not valid.'], scanner_parameters.curr_dataset);
    velocity = [];
    number_frames = [];
    return;
end

curr_colorbarbounds = scanner_parameters.curr_colorbarbounds; % Extract current colorbar bounds
color_bar = scanner_parameters.cmap; % Extract current colormap
dataset = scanner_parameters.curr_dataset;

m_img = img_size.m_img; % Number of rows in the image
n_img = img_size.n_img; % Number of columns in the image
number_frames = length(listing)-2; % Calculate number of frames in cine loop

velocity = zeros(m_img,n_img,number_frames); % Initialize velocity map matrix

% f = waitbar(number_frames,{sprintf('Calculating Velocity Map - %s',dataset),[sprintf('%0.1f',0),'%']}); % Initialize waitbar

for frame = 3:number_frames+2
    Vmap = velocity(:,:,frame-2);
%     disp(frame-2);
    perc = (frame-2)/number_frames; % Fraction complete
%     waitbar(perc,f,{sprintf('Calculating Velocity Map - %s',dataset),[sprintf('%0.1f',perc*100),'%']}) % Update wait bar
    
    slash = foldername(end-5);
    filename = [foldername slash listing(frame).name];
    img = imread(filename);
    
    velocity_spectrum = linspace(curr_colorbarbounds,-curr_colorbarbounds,size(color_bar,1)); %[cm/s]
    
    % this is vectorized version of original function
    img_trim = double(img(min(ROI_idx(:,2)):max(ROI_idx(:,2)), min(ROI_idx(:,1)):max(ROI_idx(:,1)),:));

    diff12 = abs(img_trim(:,:,1)-img_trim(:,:,2));
    diff23 = abs(img_trim(:,:,2)-img_trim(:,:,3));
    diff13 = abs(img_trim(:,:,1)-img_trim(:,:,3));

    sum_diff = diff12+diff23+diff13;    
    
    grayscale = ((img_trim(:,:,1) == img_trim(:,:,2)) & (img_trim(:,:,2) == img_trim(:,:,3)) & (img_trim(:,:,1) == img_trim(:,:,3)));
    % end of vectorization of previous


    img_color = zeros(size(img_trim));    
    
    % find all non-grayscale (and other filtered) linear (in 1 color) indices/pixels
%     inds = find(~grayscale & sum_diff>40 & diff12<20 & ~(diff23<60 & abs(img_trim(:,:,1))>40));
    inds = find(~grayscale & sum_diff>40 & diff12<20);

    % linear indices in 3 color as well
    color_inds = [inds; inds+numel(grayscale); inds+2*numel(grayscale)];
    img_color(color_inds) = img_trim(color_inds);
    
    V_trim = zeros(size(grayscale));

    
    if strcmp(scanner_parameters.curr_img_mode, 'SMI')
        % real piecewise, uses thresh, next-intercept, and coefficients
        % manually copied from linfitting script. Could automate but meh
        thresh = [12,11.444444444444445,40.333333333333336];
        
        % define images for each part of piecewise function.
        l = img_color(:,:,1)<thresh(1) & img_color(:,:,2)<thresh(2) & img_color(:,:,3)<thresh(3) & img_color(:,:,3)>0;
    %     m = img_color(:,:,3)<h_thresh & img_color(:,:,3)>=l_thresh;   % not currently using middle inds, could be useful for other modes
        h = img_color(:,:,1)>=thresh(1) | img_color(:,:,2)>=thresh(2) | img_color(:,:,3)>=thresh(3);
    
        % find pixel (linear 1 color) indices for each part of the piecewise function
        lowerinds = find(l);
    %     midinds = find(m);    % not currently using middle inds, could be useful for other modes
        upperinds = find(h);
        
    
        % calculate v from linear models (note conversion of linear 1 color to
        % linear 3 color indices)
        V_trim(lowerinds) = 0.00066101*img_color(lowerinds) + ...
                0.0032437*img_color(lowerinds+numel(grayscale)) + 0.0037869*img_color(lowerinds+2*numel(grayscale));
    
    %     not currently using middle inds, could be useful for other modes
    %     V_trim(midinds) = 0.16696  + 0.00089488*img_color(midinds+2*numel(grayscale));
    
        V_trim(upperinds) = 0.1978 + 0.00068578*img_color(upperinds) + ...
                0.0036645*img_color(upperinds+numel(grayscale)) + 1.3577e-05*img_color(upperinds+2*numel(grayscale));
    
        % clips v at max velocity
        V_trim = min(V_trim, 1);
    
    % CHI function
    else
        V_trim = -1.016 + 0.00029411*img_trim(:,:,1) + ...
                0.003738*img_trim(:,:,2) + 0.0037767*img_trim(:,:,3);
        V_trim = min(V_trim, 1);
        V_trim = max(V_trim, -1);
        V_trim = (V_trim+1)/2;
    end

    % plots for testing
%     subplot(3,1,1), plot(V_trim(inds), img_color(inds), 'LineStyle', 'none', 'Marker','.');
%     subplot(3,1,2), plot(V_trim(inds), img_color(inds+numel(grayscale)), 'LineStyle', 'none', 'Marker','.');
%     subplot(3,1,3), plot(V_trim(inds), img_color(inds+2*numel(grayscale)), 'LineStyle', 'none', 'Marker','.');


    % convert to actual max velocity
    velocity(min(ROI_idx(:,2)):max(ROI_idx(:,2)), min(ROI_idx(:,1)):max(ROI_idx(:,1)),frame-2) = V_trim*curr_colorbarbounds;
end
% close(f); % Close waitbar
end