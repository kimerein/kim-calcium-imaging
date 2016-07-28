% Script reads in output of CNMF and analyzes effects of opto stim on Ca2+
% traces

%% Set locations to files and directories
orchestraOutput='\\research.files.med.harvard.edu\Neurobio\MICROSCOPE\Kim\Data from Imaging Rig\Sabatini ScanImage Data\20160610\mouse 1\Best CNMF output\Components 200\cnmf.out';
acq_obj='\\research.files.med.harvard.edu\Neurobio\MICROSCOPE\Kim\Data from Imaging Rig\Sabatini ScanImage Data\20160610\mouse 1\a2astim.mat';
saveDir='\\research.files.med.harvard.edu\Neurobio\MICROSCOPE\Kim\Data from Imaging Rig\Sabatini ScanImage Data\20160610\mouse 1\Best CNMF output\Components 200';

%% Read in CNMF output
[Yr,b2,f2,Cn,Yk,Cf,Df,Ao]=readOrchestraOutput(orchestraOutput);

%% Fill in options used for CNMF
K = 200;                                           % number of components to be found
tau = 5;                                          % std of gaussian kernel (size of neuron)
p = 1;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.75;                                  % merging threshold
options = CNMFSetParms(...
'd1',128,'d2',512,...                         % dimensions of datasets
'ssub',6,'tsub',1,...
'min_corr',0.3,...
'nb',1,...
'save_memory',1,...
'init_method','greedy',...
'search_method','ellipse','dist',3,...      % search locations when updating spatial components
'deconv_method','constrained_foopsi',...    % activity deconvolution method
'temporal_iter',2,...                       % number of block-coordinate descent steps
'fudge_factor',0.98,...                     % bias correction for AR coefficients
'merge_thr',merge_thr,...                   % merging threshold
'gSig',tau...
);

%% Load Acquisition2P object (output of motion correction)
a=load(acq_obj);
names=fieldnames(a);
obj=a.(names{1});

%% Run analysis Step 1
dFoverF_viewer(Cf,obj,'Opto_Stim',Yk,Ao,Cn,b2,f2,Df,options,'Wheel_Encoder',saveDir,[]);

%% Run analysis Step 2
withinCellResponses=analysisSecondHalf([saveDir '\partwayData_moviematched']);

%% Run analysis Step 3