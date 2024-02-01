%%

d = string({listing.bounds});
d = [d "All match!"];
[indx,tf] = listdlg('PromptString',{'Select any values that don''t match the  image.',...
    'Empty rows correspond to 0.',''},...
    'SelectionMode','multiple','ListString',d);

if length(indx) > 1 || ~strcmp(d(indx(1)), d(end))
    prompt = {};

    for ii = 1:length(indx)
        prompt(end+1) = {['Enter correct value for row ', num2str(indx(ii)), ...
                          '. Incorrect value: ', char(d(indx(ii)))]};
    end

    corrs = str2num(char(inputdlg(prompt)));

    for ii = 1:length(indx)
        listing(indx(ii)).bounds = corrs(ii);
    end
end

%%

figure; imshow(listing(29).image);

objectRegion=round(getPosition(imrect));
close;

objectRegion;


%% Get info for selected DICOM

info = dicominfo(fullfile(fpath, listing(25).name));


%% Play selected movie

[X, cmap, alpha, overlays] = dicomread(fullfile(fpath, listing(3).name));
% implay(X);

%% Display all 1st frames and selected 1st frame

close all

I = listing(11).image;
height = length(I(:,1,1));
width = length(I(1,:,1));

imshow(I);


figure
set(gcf, 'Position',  [1600, 100, width, height])
for ii = 1:length(listing)
    
    subplot(ceil(length(listing)/6), 6, ii), imshow(listing(ii).image);
end


%% Add Mode field to pig_test_dicom_info_2.mat


modes = {'BMODE'; 'CD'; 'PD'; 'ADF'; 'ADF'; 'SMI'; 'BMODE'; 'CD'; 'PD'; 'ADF'; 'SMI'; 'BMODE'; 'CD'; 'PD'; 'ADF'; 'SMI'; 'BMODE'; 'CD'; 'PD'; 'ADF'; 'SMI'};

for ii =1:length(listing)
    listing(ii).ModeCheck = modes(ii);
end

%%

% THE FOLLOWING SECTIONS HAVE BEEN ADDED INTO FUNCTIONS export_dicom_frames
% AND extract_dicom_info. USE OF THOSE IS PREFERRED

%% Load mat 

fpath = 'C:\Users\Denis\Documents\JHSOM\Long-Term Storage PhD\Data\220111 SCAR initial\2022_1_11_Ultrasound';
% fpath = 'C:\Users\Denis\OneDrive for Business\HEPIUS\Pig\Experiments\Terminal_Ultrasound_2021_12_16';
% fpath = 'C:\Users\Denis\Documents\JHSOM\Long-Term Storage PhD\Data\Test Data for MATLAB';

fpath = uigetdir;

fnames = dir(fpath);
inds = [];

warning off images:dicomparse:shortImport

f = waitbar(0, 'Extracting DICOM File Data...');

for ii = 1:length(fnames)
        
    if ~contains(fnames(ii).name(1), 'A')
        continue
    end
    
    waitbar(ii/length(fnames),f,['Loading DICOM Files (' fnames(ii).name ')...']);

    inds = [inds ii];

    info = dicominfo(fullfile(fpath, fnames(ii).name));
    
    fnames(ii).AcquisitionDate = info.AcquisitionDate;

    fnames(ii).timestring = info.AcquisitionTime;
    fnames(ii).timestamp = 3600*str2double(info.AcquisitionTime(1:end-8)) ...
        + 60*str2double(info.AcquisitionTime(end-7:end-6))...
        + str2double(info.AcquisitionTime(end-5:end));

    try
        fnames(ii).NumberOfFrames = info.NumberOfFrames;
        fnames(ii).FrameTime = info.FrameTime;
        fnames(ii).EffectiveDuration = info.EffectiveDuration;
        fnames(ii).CineRate = info.CineRate;
%         fnames(ii).ImageType = info.ImageType(end-3:end);
%         fnames(ii).testvar = info.HighBit;
    catch
        fnames(ii).NumberOfFrames = 1;
        fnames(ii).FrameTime = 0;
        fnames(ii).EffectiveDuration = 0;
        fnames(ii).CineRate = 0;
%         fnames(ii).ImageType = 0;
%         fnames(ii).testvar = 0;
    end
    
    fnames(ii).image = dicomread(info,"Frames",1);
    
end

warning on images:dicomparse:shortImport

close(f);

listing = fnames(inds);
height = length(listing(1).image(:,1,1));
width = length(listing(1).image(1,:,1));


%% Setting the clock to 0 at first DICOM

t_0 = listing(1).timestamp;

for ii = 1:length(listing)
    listing(ii).time = listing(ii).timestamp - t_0;    
end


%% Recognize Modality and colorbar bounds

close all;

figure
set(gcf, 'Position',  [1600, 100, width, height])


for ii = 1:length(listing)
    
    graphRegion = imresize(listing(ii).image(775:790, 1:27, :), 2);
    graphRegion = imbinarize(rgb2gray(graphRegion));
    axisLabel = ocr(graphRegion, 'TextLayout','Word', 'CharacterSet','cm');
    
    if contains(axisLabel.Text, 'cm')
        listing(ii).Mode = "spectral";
        warning([listing(ii).name ' is a spectral image/video. It will not work for spatial flow analysis.'])
        continue;
    end

    colorlist = {'ADF', 'ADF', 'CDI', 'CDI', 'PDI', 'SMI', 'BMODE', 'SWE'};
    colorvalues = [240 220 220; 220 220 240; 0 190 255; 255 190 0; 252 224 0; 162 161 254; 230 230 230; 175, 0, 0];
    
    colorsample = double(listing(ii).image(186,1218,:));

    error = (colorvalues(:,1) - colorsample(1)).^2 + (colorvalues(:,2) - colorsample(2)).^2 + ...
            (colorvalues(:,3) - colorsample(3)).^2;
    [M, colorind] = min(error);
    
    listing(ii).error = M;

    listing(ii).Mode = string(colorlist(colorind));


    if listing(ii).Mode == "BMODE"
        listing(ii).bounds = [];
    
    else
        boundsRegion = imresize(listing(ii).image(155:170, 1195:1225, :), 2);
    
        boundsRegion = imbinarize(rgb2gray(boundsRegion));
    
        subplot(ceil(length(listing)/6), 6, ii), imshow(boundsRegion);
        cbarBounds = ocr(boundsRegion, 'TextLayout','Word', 'CharacterSet','.0123456789');
    
        listing(ii).bounds = cbarBounds.Text;

        if listing(ii).Mode == "SWE" && isempty(listing(ii).bounds)
            boundsRegion = imresize(listing(ii).image(155:170, 42:72, :), 2);
            boundsRegion = imbinarize(rgb2gray(boundsRegion));

            subplot(ceil(length(listing)/6), 6, ii), imshow(boundsRegion);
            cbarBounds = ocr(boundsRegion, 'TextLayout','Word', 'CharacterSet','.0123456789');
        
            listing(ii).bounds = cbarBounds.Text;
        end
    end

end


%% Export the nth (iith) frame of each DICOM (to run in Kelley's program)

newfolder = ['Frame ' num2str(ii)];
date = listing(1).AcquisitionDate;

if ~exist(fullfile(fpath, date,newfolder),'dir')
    mkdir(fullfile(fpath, date,newfolder));
end

warning off images:dicomparse:shortImport
for dicomTag = 1:length(listing)
    fname = fullfile(fpath, listing(dicomTag).name);

    number_frames = listing(dicomTag).NumberOfFrames;
    dataset = listing(dicomTag).name;

    ii = 1;
    info = dicominfo(fullfile(fpath, listing(dicomTag).name));
    image = dicomread_debugged(info,"Frames",ii);

    if ~exist(fullfile(fpath, date,newfolder,[dataset '.png']),'dir')
        disp('Saving image');
        imwrite(image, fullfile(fpath, date,newfolder,[dataset '.png']),'png');
%         save(fullfile(fpath, date,'Images',dataset,['image' num2str(ii)]), image);

    else
        disp('image already exists');
    end
end
warning off images:dicomparse:shortImport

%% Export every nth frame to correct folder (to run in Kelley's program)

% Currently loks like 1500 frames will take up 5.6ish GB, which is my max
% this is 25 sec at 60 fps and 42 sec at 35 fps


dicomTag = 3;
interval = 2;

number_frames = listing(dicomTag).NumberOfFrames;
date = listing(dicomTag).AcquisitionDate;
dataset = listing(dicomTag).name;

fname = fullfile(fpath, listing(dicomTag).name);

if ~exist(fullfile(fpath, date,'Images',dataset),'dir')
    mkdir(fullfile(fpath, date,'Images',dataset));
end

warning off images:dicomparse:shortImport
for ii = 1:3 %interval:number_frames
    
    info = dicominfo(fullfile(fpath, listing(dicomTag).name));
    image = dicomread_debugged(info,"Frames",ii);

    if ~exist(fullfile(fpath, date,'Images',dataset,['image' num2str(ii) '.png']),'dir')
        disp('Saving image');
        imwrite(image, fullfile(fpath, date,'Images',dataset,['image' num2str(ii) '.png']),'png');
%         save(fullfile(fpath, date,'Images',dataset,['image' num2str(ii)]), image);

    else
        disp('image already exists');
    end
end
warning off images:dicomparse:shortImport