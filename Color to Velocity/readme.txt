Written: 1/26/2022 (Denis Routkevitch, droutke1@jhmi.edu)


NOTE: this currently only works for Canon exported DICOMs

Basic instructions:

1. Put all DICOM files in one folder.
2. Run prepreprocessing script, navigate to your DICOMs folder and select any of the files
	-Even though you select one file, all DICOMs will be analyzed.
	-To select certain DICOMs for analysis, put them in their own folder
3. Confirm that the numbers are read correctly, contact me (email at top) if not.
	- Only if troubleshooting is enable in the prepreprocessing script
4. Run preprocessing script.
	-Navigate to DICOMs folder, then select "DICOM listing.mat"
	-Select the ROIs.
5. Run the main script
	-Navigate to DICOMs folder, then to \**date**\PreProcessing, and select any of the *.mat files
6. Run the postprocessing script
	-Navigate to DICOMs folder, then to \**date**\Velocity Maps, and select any of the *.mat files

Velocity maps and quantitative data folders contain *.mat files ready for further analysis!

Don't hesitate to email me with further questions!