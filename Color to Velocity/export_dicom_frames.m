% Function: export_dicom_frames
%
% Purpose: Export all images
%
% Input parameters:
%   listing: struct (same as full_dicom_listing)
%   interval: double (whole number, downsample frames by this amount)
%   dicom_inds: double or double array, choose which dicoms to export
%           Can also accept char 'all' to export all dicoms
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
%
%
% Created by: Denis Routkevitch (droutke1@jhmi.edu)
% Edited: 1/28/2022
%       Increased speed of export via grouped loading from DICOM 
% Edit 2/5/2022  (DR):
%       Changed naming system
% Edit 4/5/2022  (DR):
%       Changed description


function export_dicom_frames(listing, interval, dicom_inds)

msg = MException('funcInput:wrongVariableType', ['export_dicom_frames(listing, interval, dicom_).' ...
        ' If included in argslist, listing must be a structure, interval must be a whole number double ' ...
        'and dicom_inds must be ''all'' or an array of whole number doubles']);

% check inputs
if nargin == 1 & isa(listing, 'struct')
    interval = 1;
    dicom_inds = 1:length(listing);
elseif nargin == 2 & isa(listing, 'struct') & isa(interval, 'double')
    dicom_inds = 1:length(listing);
elseif nargin == 3 & isa(listing, 'struct') & isa(interval, 'double') & isa(dicom_inds, 'double')
elseif nargin == 3 & isa(listing, 'struct') & isa(interval, 'double') & isa(dicom_inds, 'char')
    if strcmpi(dicom_inds, 'all')
        dicom_inds = 1:length(listing);
    else
        throw(msg);
    end
else    
    throw(msg);
end

fpath = listing(1).folder;
save(fullfile(fpath, 'DICOM listing.mat'), 'listing');


f = waitbar(0, 'Exporting DICOM Files...');    % initialize waitbar
for dicomTag = dicom_inds
    number_frames = listing(dicomTag).NumberOfFrames;
    date = listing(dicomTag).AcquisitionDate;
    dataset = listing(dicomTag).name;
    
    fpath = listing(dicomTag).folder;
    fname = fullfile(fpath, listing(dicomTag).name);
    
    % check if Images folder exists
    if ~exist(fullfile(fpath, date,'Images',dataset),'dir')
        mkdir(fullfile(fpath, date,'Images',dataset));
    end

    % initialize waitbar
    k = find(dicom_inds==dicomTag, 1);
    waitbar(0,f,sprintf('Exporting DICOM File %s (%f/%f)', listing(dicomTag).name,...
        k,length(dicom_inds)));
    
    warning off images:dicomparse:shortImport
    batch = 10;
    % Load DICOM frames in batches
    for ii = 1:interval*batch:number_frames
        
        waitbar(ii/number_frames,f,sprintf('Exporting DICOM File %s (%g/%g)', listing(dicomTag).name,...
            k,length(dicom_inds)));
        
        frameset  = ii:interval:min(ii+interval*batch-1, number_frames);
        
        % check for existance of first frame in batch. Lazy solution, coulf
        % be difficulties if run multiple times with different batch sizes
        if ~exist(fullfile(fpath, date,'Images',dataset,sprintf('image%04d.png', ii)),'file')
            
            % read selected frames
            info = dicominfo(fullfile(fpath, listing(dicomTag).name));
            image = dicomread_debugged(info,"Frames",frameset);
            
            for jj = 1:length(frameset)
                % export image as png
                imwrite(image(:,:,:,jj), fullfile(fpath, date,'Images',dataset,sprintf('image%04d.png', frameset(jj))),'png');
        %         save(fullfile(fpath, date,'Images',dataset,['image' num2str(ii)]), image);
                
                
            end
        end
    end
    
end
warning on images:dicomparse:shortImport

close(f);

end