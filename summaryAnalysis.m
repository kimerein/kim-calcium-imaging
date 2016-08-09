function summaryAnalysis(dataDir,saveDir,acq_obj_pointer)

% Set locations to files and directories
orchestraOutput=dataDir;
acq_obj=acq_obj_pointer;
saveDir=saveDir;

% Read in CNMF output
[Yr,b2,f2,Cn,Yk,Cf,Df,Ao]=readOrchestraOutput(orchestraOutput);

% Fill in options used for CNMF
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

% Load Acquisition2P object (output of motion correction)
a=load(acq_obj);
names=fieldnames(a);
obj=a.(names{1});



% Iterate through analysis conditions

% Clear counter
iterateAnalysisSettings([],[],[],[],[],[],[],1,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in saved data
load([saveDir '\partwayData_moviematched\optoMapping.mat']);
load([saveDir '\partwayData_moviematched\optoStimTypes.mat']);
load([saveDir '\partwayData_moviematched\useComponents.mat']);
readin_data.optoMapping=optoMapping;
readin_data.optoStimTypes=optoStimTypes;
readin_data.useComponents=useComponents;

% First, real opto stim
dFoverF_viewer(Cf,obj,'Opto_Stim',Yk,Ao,Cn,b2,f2,Df,options,'Wheel_Encoder',saveDir,readin_data);

% % Load in saved data
% load([saveDir '\partwayData_moviematched\optoMapping.mat']);
% load([saveDir '\partwayData_moviematched\optoStimTypes.mat']);
% load([saveDir '\partwayData_moviematched\useComponents.mat']);
% readin_data.optoMapping=optoMapping;
% readin_data.optoStimTypes=optoStimTypes;
% readin_data.useComponents=useComponents;

% Run analysis on real opto stim
[withinCellResponses,withinCellStats,withinCellAverages,times,optoForProfile,matrixOfEffects]=analysisSecondHalf([saveDir '\partwayData_moviematched']);
traceOutput=plotCaResponse(withinCellAverages,withinCellStats,times,optoForProfile);
fxDistributionOutput=analyzeTrialByTrial(withinCellResponses,withinCellStats,times);

% Save output of analysis
if ~exist([saveDir '\opto_stim'],'dir')
    mkdir([saveDir '\opto_stim']);
end
save([saveDir '\opto_stim\withinCellResponses.mat'],'withinCellResponses');
save([saveDir '\opto_stim\withinCellStats.mat'],'withinCellStats');
save([saveDir '\opto_stim\withinCellAverages.mat'],'withinCellAverages');
save([saveDir '\opto_stim\times.mat'],'times');
save([saveDir '\opto_stim\optoForProfile.mat'],'optoForProfile');
save([saveDir '\opto_stim\matrixOfEffects.mat'],'matrixOfEffects');
save([saveDir '\opto_stim\traceOutput.mat'],'traceOutput');
save([saveDir '\opto_stim\fxDistributionOutput.mat'],'fxDistributionOutput');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change analysis settings to shutter only (no opto stim)
iterateAnalysisSettings([],[],[],[],[],[],[],0,1);

% Shutter only
dFoverF_viewer(Cf,obj,'Opto_Stim',Yk,Ao,Cn,b2,f2,Df,options,'Wheel_Encoder',saveDir,readin_data);

% Run analysis on shutter only
[withinCellResponses,withinCellStats,withinCellAverages,times,optoForProfile,matrixOfEffects]=analysisSecondHalf([saveDir '\partwayData_moviematched']);
traceOutput=plotCaResponse(withinCellAverages,withinCellStats,times,optoForProfile);
fxDistributionOutput=analyzeTrialByTrial(withinCellResponses,withinCellStats,times);

% Save output of analysis
if ~exist([saveDir '\control'],'dir')
    mkdir([saveDir '\control']);
end
save([saveDir '\control\withinCellResponses.mat'],'withinCellResponses');
save([saveDir '\control\withinCellStats.mat'],'withinCellStats');
save([saveDir '\control\withinCellAverages.mat'],'withinCellAverages');
save([saveDir '\control\times.mat'],'times');
save([saveDir '\control\optoForProfile.mat'],'optoForProfile');
save([saveDir '\control\matrixOfEffects.mat'],'matrixOfEffects');
save([saveDir '\control\traceOutput.mat'],'traceOutput');
save([saveDir '\control\fxDistributionOutput.mat'],'fxDistributionOutput');
