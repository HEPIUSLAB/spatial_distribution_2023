% Function: extract_dicom_info
%
% Purpose: Open DICOM file metadata and extract relevant information
% Then, parse image for info on colorbar bounds and Doppler modality
%
% Input parameters:
%   File path (optional): char (default is uigetdir)
%   Check colorbars (optional): logical (default false)
%
% Output parameters:
%   listing: struct (directory information+metadata)
%       listing.name: char
%       listing.folder: char
%       listing.date: char
%       listing.bytes: double
%       listing.isdir: logical
%       listing.datenum: double
%       listing.AcquisitionDate: char
%       listing.time: double (in seconds with 0 being start of first
%           recording
%       listing.NumberOfFrames: double
%       listing.FrameTime: double (in milliseconds)
%       listing.image: uint8 (first frame of DICOM)
%       listing.Mode: string
%       listing.ColorBarBounds: double
%       listing.delxy: double
%
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)
% Edited: 2/2/2022
%       Solved indexing bug and added poor colorbar bounds recognition detection
%       Resolved issue of dialog box too large in some cases
%       Resolved issue of double digit colorbar bounds
% Last edited: 3/2/2022
%       added delxy as parameter to extract
% Edit 4/5/2022 (DR):
%       Changed frame export method (last frame)

function listing = extract_dicom_info(fpath, checkCBar)

% check inputs
if nargin == 0
    [~, fpath] = uigetfile('*.*');
    checkCBar = 0;
elseif nargin == 1 & isa(fpath, 'char')
    checkCBar = 0;
elseif nargin == 1 & isa(fpath, 'logical')
    checkCBar = fpath;
    [~, fpath] = uigetfile('*.*');
elseif nargin == 2 & isa(fpath, 'char') & isa(checkCBar, 'logical')
else
    msg = MException('funcInput:wrongVariableType', ['extract_dicom_info(fpath, checkCBar). If included in argslist,' ...
        'fpath (file path) must be a char and checkCBar must be logical (true or false)']);
    throw(msg);
end

fnames = dir(fpath);    % initialize filenames
inds = [];


warning off images:dicomparse:shortImport   % this warning is annoying and harmless

f = waitbar(0, 'Extracting DICOM File Data...');    % initialize waitbar

for ii = 1:length(fnames)
        
    if ~contains(fnames(ii).name(1), 'A')   % only work on A____ files
        continue
    end
    
    % update waitbar
    waitbar(ii/length(fnames),f,['Loading DICOM Files (' fnames(ii).name ')...']);

    inds = [inds ii];   % indices of DICOMs

    info = dicominfo(fullfile(fpath, fnames(ii).name)); % load metadata
    
    fnames(ii).AcquisitionDate = info.AcquisitionDate;  % get surgery date

    fnames(ii).timestring = info.AcquisitionTime;   % raw timestamp

    
    
    % convert raw time to seconds
    fnames(ii).timestamp = 3600*str2double(info.AcquisitionTime(1:end-8)) ...
        + 60*str2double(info.AcquisitionTime(end-7:end-6))...
        + str2double(info.AcquisitionTime(end-5:end));
    
    % extract movie specific parameters and set to defaults if image
    try
        fnames(ii).NumberOfFrames = info.NumberOfFrames;
        fnames(ii).FrameTime = info.FrameTime;
        fnames(ii).delxy = info.SequenceOfUltrasoundRegions.Item_1.PhysicalDeltaX;  % cm/pixel
%         fnames(ii).EffectiveDuration = info.EffectiveDuration;
%         fnames(ii).CineRate = info.CineRate;
%         fnames(ii).ImageType = info.ImageType(end-3:end);
%         fnames(ii).testvar = info.HighBit;
    catch
        fnames(ii).NumberOfFrames = 1;
        fnames(ii).FrameTime = 0;
        fnames(ii).delxy = [];  % cm/pixel
%         disp(fnames(ii).name)
%         fnames(ii).EffectiveDuration = 0;
%         fnames(ii).CineRate = 0;
%         fnames(ii).ImageType = 0;
%         fnames(ii).testvar = 0;
    end
    
    % export image frame (last if possible)
    try
        fnames(ii).image = dicomread_debugged(info,"Frames",info.NumberOfFrames);  % Extract last frame
    catch
        fnames(ii).image = dicomread_debugged(info,"Frames",1);     % extract first frame if that fails
    end
end

% turn warning back on
warning on images:dicomparse:shortImport


% remove non-DICOM files from listing
listing = fnames(inds);
% height = length(listing(1).image(:,1,1));
% width = length(listing(1).image(1,:,1));


%%%% set time of first recording as t = 0 %%%%
t_0 = listing(1).timestamp;

for ii = 1:length(listing)
    listing(ii).time = listing(ii).timestamp - t_0;    
end


% Remove unnecessary fields
listing = rmfield(listing, 'timestring');
listing = rmfield(listing, 'timestamp');
listing = rmfield(listing, 'datenum');
listing = rmfield(listing, 'isdir');


%%%% Recognize Modality and colorbar bounds %%%%
close all;


% initialize figure for loading colorbar bound images

rows = 15;
cols = ceil(length(listing)/rows);
plotind = reshape(1:(rows*cols), cols, rows).';

if checkCBar
    f2 = figure;
    screen = get(0,'ScreenSize');
    figheight = min([screen(4)-150 750]);
    set(gcf, 'Position',  [200, 50, 100*cols, figheight])
end


for ii = 1:length(listing)
    
    waitbar(ii/length(listing),f,['Extracting DICOM Info (' listing(ii).name ')...']); % update waitbar
    
    % check for spectral analysis (not suitable for further input)
    % isolate image at coordinates where 'cm' appears
    graphRegion = imresize(listing(ii).image(775:790, 1:27, :), 2);
    graphRegion = imbinarize(rgb2gray(graphRegion));
    axisLabel = ocr(graphRegion, 'TextLayout','Word', 'CharacterSet','cm');
    
    % if there is a "cm" axis label, image is spectral
    if contains(axisLabel.Text, 'cm')
        listing(ii).Mode = "spectral";
        listing(ii).bounds = 0; % not calculating spectral bounds for now
        
        % show blank image if enabled in argin
        if checkCBar
            subplot(rows, cols, plotind(ii)), imshow((zeros(320, 620)));
            ylabel(listing(ii).name, 'Rotation', 0, 'FontSize', 6);
        end

        warning([listing(ii).name ' is a spectral image/video. It will not work for spatial flow analysis.'])
        continue;
    end


    % Colorbar values for each modality at the extraction coordinate (186, 1218)
    colorlist = {'ADF', 'ADF', 'CDI', 'CDI', 'PDI', 'SMI', 'BMODE', 'SWE', 'CHI'};
    colorvalues = [240 220 220; 220 220 240; 0 190 255; 255 190 0; 252 224 0; 162 161 254; 230 230 230; 175, 0, 0; 254 246 233];
    
    % extract color bar color and calculate sum of squared error for each
    % modality
    colorsample = double(listing(ii).image(186,1218,:));
    error = (colorvalues(:,1) - colorsample(1)).^2 + (colorvalues(:,2) - colorsample(2)).^2 + ...
            (colorvalues(:,3) - colorsample(3)).^2;     % calculate sum of squared error
    
    % minimize error to determine modality
    [~, colorind] = min(error);
    listing(ii).Mode = string(colorlist(colorind));

    
    % no bounds if B-mode, fixed, non symmetric bounds for SWE, handle
    % later if presets for SWE change

    if listing(ii).Mode == "BMODE" || listing(ii).Mode == "SWE" || listing(ii).Mode == 'CHI'
        listing(ii).bounds = 0;
        
        % show blank image if enabled in argin
        if checkCBar
            subplot(rows, cols, plotind(ii)), imshow((zeros(320, 620)));
            ylabel(listing(ii).name, 'Rotation', 0, 'FontSize', 6);
        end
    else
        % extract and process region where colorbar bounds appear
        boundsColorRegion = imresize(listing(ii).image(155:170, 1188:1225, :), 2);    
        boundsRegion = imbinarize(rgb2gray(boundsColorRegion));
        
        
        % show bounds image if enabled in argin
        if checkCBar
            subplot(rows, cols, plotind(ii)), imshow(boundsColorRegion);
            ylabel(listing(ii).name, 'Rotation', 0, 'FontSize', 6);
        end
        
        % read cbar value and enter into listing
        cbarBounds = ocr(boundsRegion, 'TextLayout','Word', 'CharacterSet','.0123456789');    
        listing(ii).bounds = str2double(cbarBounds.Text);
        
        % check for poor reading
        if isnan(listing(ii).bounds) || min(cbarBounds.CharacterConfidences) < 0.8
            listing(ii).bounds = NaN;
            warning([listing(ii).name ' has poor number recognition. Check bounds for accuracy!!'])
        end
    end

end

close(f);   % close waitbar

% check cbar bound values if enabled
if checkCBar

    % create options and open dialog box
    d = string({listing.bounds});
    
    for ii = 1:length(d)
        if ismissing(d(ii))
            d(ii) = "No value detected";
        end
        
        d(ii) = strcat(d(ii), " (", string(listing(ii).name), ")");
    end

    d = ["All match!" d];
    [indx,tf] = listdlg('PromptString',{'Select any values that don''t match the  image.',...
        'Empty rows correspond to 0.',''},...
        'SelectionMode','multiple','ListString',d);
    
    % if dialog box is anything but 'All match!', prompt for corrections
    if length(indx) > 1 || indx(1)~=1
        for kk = 1:5:max(length(indx)-1,1)
            prompt = {};
            
            if indx(1) == 1
                indx = indx(2:end);
            end
    
            for ii = indx(kk:min(kk+4,length(indx)))
                prompt(end+1) = {sprintf(['Enter correct value for row ', num2str(ii-1), ...
                                  ' (%s). Incorrect value: ', char(d(ii))], listing(ii-1).name)};
            end
        
            corrs = str2num(char(inputdlg(prompt)));
            
            % implement corrections
            for jj = kk:min(kk+4,length(indx))
                listing(indx(jj)-1).bounds = corrs(jj-kk+1);
            end
        end
    end
    
    close(f2); % close images of bounds
end

end