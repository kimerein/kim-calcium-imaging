function extractSingleTrialData(datadir)

datadir='\\research.files.med.harvard.edu\neurobio\MICROSCOPE\Kim\Data from Imaging Rig\Sabatini ScanImage Data\Results\';

load([datadir 'dir.mat']);
load([datadir 'optoStimTypes_forExpt.mat']);
load([datadir 'allexpts_withTrialByTrialBeh.mat']);
load([datadir 'opto_aligned_data.mat']);
load([datadir 'useTheseOpto_forDifferentOptoTypes.mat']);

[matrixCellResponses,matrixBeh,selectCellResponses,selectBeh,allBehProfiles,allTrialsResponses,allTrialsBeh,alloptoprofiles]=selectDataByBehavior(allSingleTrialBeh,allTrialResponses{24},allOptoData{24},timestep,allOptoTypes,useTheseOpto,allexpts.cell_from_expt{24});

[fromOpto_cellResponses,fromOpto_beh]=selectFinalDatafromBehAndOpto(allTrialsResponses,allTrialsBeh,alloptoprofiles,timestep,allOptoData{24},opto20Hz_4V,allexpts.cell_from_expt{24});

si=plotZscoredData(fromOpto_cellResponses,0:timestep:(size(allSingleTrialBeh{1},2)-1)*timestep,'opto20Hz_4V',fromOpto_beh,[]);