function plot_raw_vs_processed(dataDir,saveDir,acq_obj_pointer,cell_by_traces)

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

if ~exist([saveDir '\comparing_raw_to_processed'],'dir')
    mkdir([saveDir '\comparing_raw_to_processed']);
end

% First, real opto stim
% cell_by_traces=(cell_by_traces*max(max(Cf)))./max(max(cell_by_traces));
% dFoverF_viewer(cell_by_traces,obj,'Opto_Stim',Yk,Ao,Cn,b2,f2,Df,options,'Wheel_Encoder',[saveDir '\comparing_raw_to_processed'],readin_data);
cell_by_traces=(cell_by_traces*max(max(Yk)))./max(max(cell_by_traces));
dFoverF_viewer(Cf,obj,'Opto_Stim',cell_by_traces,Ao,Cn,b2,f2,Df,options,'Wheel_Encoder',[saveDir '\comparing_raw_to_processed'],readin_data);

end