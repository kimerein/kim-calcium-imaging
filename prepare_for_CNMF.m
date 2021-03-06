% Script prepares raw data for CNMF analysis

%% Set variables
saveDir='\\research.files.med.harvard.edu\Neurobio\MICROSCOPE\Kim\Data from Imaging Rig\Sabatini ScanImage Data\20160607';

%% Make Acquisition2P object

obj=Acquisition2P([],@SC2Pinit);

%% Motion correction

obj.motionCorrect;
save([saveDir '\obj.mat'],'obj');

%% Combine tif movies

combineAndDivideTifs(1,50,12);

%% Run CNMF on Orchestra