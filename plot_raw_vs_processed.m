function plot_raw_vs_processed(dataDir,saveDir,acq_obj_pointer,cell_by_traces)

[Selection,ok]=listdlg('PromptString','Select a comparison',...
                       'SelectionMode','single',...
                       'ListString',{'raw vs backgnd-subtracted','raw vs CNMF inferred'});
if ok==0
    return
end

% Set locations to files and directories
orchestraOutput=dataDir;
acq_obj=acq_obj_pointer;

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

if ~exist([saveDir '\comparing_raw_to_processed'],'dir')
    mkdir([saveDir '\comparing_raw_to_processed']);
end

if isempty(cell_by_traces)
    % Get ROIs from CNMF
    roi_masks=get_ROIs(Ao,options);
    
    % Normalize ROIs such that spatial integral across all pixels is 1
    for i=1:length(roi_masks)
        temp=roi_masks{i};
        integral=sum(sum(temp));
        roi_masks{i}=temp/integral;
    end
    
    % Get Ca2+ traces by applying ROI masks to movie
    traces=apply_ROIs_to_movie(roi_masks,acq_obj);
    
    % Get traces concatenated across trials
    cell_by_traces=convert_trials_to_continuous(traces);
    
    save([saveDir '\comparing_raw_to_processed\cell_by_traces.mat'],'cell_by_traces');
end

% Load Acquisition2P object (output of motion correction)
a=load(acq_obj);
names=fieldnames(a);
obj=a.(names{1});

% Load in saved data
load([saveDir '\partwayData_moviematched\optoMapping.mat']);
load([saveDir '\partwayData_moviematched\optoStimTypes.mat']);
load([saveDir '\partwayData_moviematched\useComponents.mat']);
readin_data.optoMapping=optoMapping;
readin_data.optoStimTypes=optoStimTypes;
readin_data.useComponents=useComponents;

% First, real opto stim
if Selection==1
    cell_by_traces=(cell_by_traces*max(max(Cf)))./max(max(cell_by_traces));
    dFoverF_viewer(cell_by_traces,obj,'Opto_Stim',Yk,Ao,Cn,b2,f2,Df,options,'Wheel_Encoder',[saveDir '\comparing_raw_to_processed'],readin_data);
elseif Selection==2
    cell_by_traces=(cell_by_traces*max(max(Yk)))./max(max(cell_by_traces));
    dFoverF_viewer(Cf,obj,'Opto_Stim',cell_by_traces,Ao,Cn,b2,f2,Df,options,'Wheel_Encoder',[saveDir '\comparing_raw_to_processed'],readin_data);
else
    return
end

[withinCellResponses,withinCellStats,withinCellAverages,times,optoForProfile,matrixOfEffects]=analysisSecondHalf([saveDir '\comparing_raw_to_processed\partwayData_moviematched']);
traceOutput=plotCaResponse(withinCellAverages,withinCellStats,times,optoForProfile);
fxDistributionOutput=analyzeTrialByTrial(withinCellResponses,withinCellStats,times);

end