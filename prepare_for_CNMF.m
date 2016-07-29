% Script prepares raw data for CNMF analysis

%% Set variables
saveDir='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\Data from Imaging Rig\Sabatini ScanImage Data\20160609';

%% Make Acquisition2P object

obj=Acquisition2P([],@SC2Pinit);
save([saveDir '\obj.mat'],'obj');

%% Motion correction

obj.motionCorrect;

%% Combine tif movies

combineAndDivideTifs(1,50,12);

%% Run CNMF on Orchestra