function roi_masks=get_CNMF_ROIs(dataDir)

% Set locations to files and directories
orchestraOutput=dataDir;

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

% Get ROIs from CNMF
roi_masks=get_ROIs(Ao,options);

% Normalize ROIs such that spatial integral across all pixels is 1
for i=1:length(roi_masks)
    temp=roi_masks{i};
    integral=sum(sum(temp));
    roi_masks{i}=temp/integral;
end

end