function export_frames_fast(listing)
fpath = listing(1).folder;
save(fullfile(fpath, 'DICOM listing.mat'), 'listing');

[~, img_all_folder] = uigetfile('*.*','Select an exported frame from Microdicom in your "Images All" folder');
export_folder = fullfile(img_all_folder,'..');%'C:\Users\kkempsk1\OneDrive - Johns Hopkins\Exporting test\Images';
listing_img = dir(img_all_folder);
date = listing(3).AcquisitionDate;

number_imgs = length(listing_img)-2;
f = waitbar(number_imgs,{'Preparing Dataset Images',[sprintf('%0.1f',0),'%']}); % Initialize waitbar

for i = 3:length(listing_img)
    perc = (i-2)/number_imgs; % Fraction complete
    waitbar(perc,f,{'Preparing Dataset Images',[sprintf('%0.1f',perc*100),'%']}) % Update wait bar

    curr_filename = listing_img(i).name;
    dataset_name = curr_filename(1:5);

    export_folder_sub = fullfile(export_folder,date,'Images',dataset_name);
    if ~exist(export_folder_sub,'dir')
        mkdir(export_folder_sub)
    end
    source_file = fullfile(img_all_folder,curr_filename);
    destination_file = fullfile(export_folder_sub,curr_filename);
    copyfile(source_file, destination_file)
end
close(f);
end