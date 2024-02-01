% Function: color_mapping_function
%
% Purpose: Create 3D matrix of the velocity map for the current dataset
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
% Last edited: 2/15/2022 (Denis Routkevitch, droutke1@jhmi.edu)
%       changed implementation slightly to permit parallel processng with
%       parfor loop.

function [velocity,number_frames] = color_mapping_function(foldername,listing,ROI_idx,img_size,scanner_parameters)

curr_colorbarbounds = scanner_parameters.curr_colorbarbounds; % Extract current colorbar bounds
color_bar = scanner_parameters.cmap; % Extract current colormap
dataset = scanner_parameters.curr_dataset;

m_img = img_size.m_img; % Number of rows in the image
n_img = img_size.n_img; % Number of columns in the image
number_frames = length(listing)-2; % Calculate number of frames in cine loop

velocity = zeros(m_img,n_img,number_frames); % Initialize velocity map matrix

f = waitbar(number_frames,{sprintf('Calculating Velocity Map - %s',dataset),[sprintf('%0.1f',0),'%']}); % Initialize waitbar

for frame = 3:number_frames+2
    Vmap = velocity(:,:,frame-2);
%     disp(frame-2);
    perc = (frame-2)/number_frames; % Fraction complete
    waitbar(perc,f,{sprintf('Calculating Velocity Map - %s',dataset),[sprintf('%0.1f',perc*100),'%']}) % Update wait bar
    
    slash = foldername(end-5);
    filename = [foldername slash listing(frame).name];
    img = imread(filename);
    
    velocity_spectrum = linspace(curr_colorbarbounds,-curr_colorbarbounds,size(color_bar,1)); %[cm/s]

    [mm,~] = size(ROI_idx);
    for m = 1:mm % Loop through each pixel in your ROI
        curr_px = squeeze(img(ROI_idx(m,2),ROI_idx(m,1),:));
        curr_px = double(curr_px);
        diff12 = abs(curr_px(1)-curr_px(2));
        diff23 = abs(curr_px(2)-curr_px(3));
        diff13 = abs(curr_px(1)-curr_px(3));
        sum_diff = diff12+diff23+diff13;
        
        grayscale = ((curr_px(1) == curr_px(2)) && (curr_px(2) == curr_px(3)) && (curr_px(1) == curr_px(3)));
        if ~grayscale
            if sum_diff > 40
                [j,~] = find((curr_px(1) == color_bar(:,:,1)) & curr_px(2) == color_bar(:,:,2) & curr_px(3) == color_bar(:,:,3));
                if ~isempty(j)
                    if numel(j)>1
                        velocity_all = velocity_spectrum(j);
                        Vmap(ROI_idx(m,2),ROI_idx(m,1)) = mean(velocity_all);
                    else
                        Vmap(ROI_idx(m,2),ROI_idx(m,1)) = velocity_spectrum(j);
                    end
                else
                    Red = squeeze(color_bar(:,:,1));
                    Green = squeeze(color_bar(:,:,2));
                    Blue = squeeze(color_bar(:,:,3));

                    R = double(curr_px(1));
                    G = double(curr_px(2));
                    B = double(curr_px(3));

                    tmp = sqrt((double(Red) - R).^2 + (double(Green) - G).^2 + (double(Blue) - B).^2);
                    [x,~] = find(tmp == min(tmp) );
                    Vmap(ROI_idx(m,2),ROI_idx(m,1)) = velocity_spectrum(mode(x));
                end
            end
        end

    end

    velocity(:,:,frame-2) = Vmap;
end
close(f); % Close waitbar
end