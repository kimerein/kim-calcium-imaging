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

function roi_masks=get_ROIs(Ao,options)

% Get options used for Ca2+ source extraction
defoptions = CNMFSetParms;
if isempty(options); options = []; end
if ~isfield(options,'d1') || isempty(options.d1); d1 = input('What is the total number of rows? \n'); else d1 = options.d1; end          % # of rows
if ~isfield(options,'d2') || isempty(options.d2); d2 = input('What is the total number of columns? \n'); else d2 = options.d2; end          % # of columns
if ~isfield(options,'plot_df') || isempty(options.plot_df); options.df = defoptions.plot_df; end
plot_df = options.plot_df;
if ~isfield(options,'make_gif') || isempty(options.make_gif); options.make_gif = defoptions.make_gif; end
make_gif = options.make_gif;
if ~isfield(options,'save_avi') || isempty(options.save_avi); options.save_avi = defoptions.save_avi; end
save_avi = options.save_avi;
if ~isfield(options,'sx') || isempty(options.sx); options.sx = defoptions.sx; end
sx = min([options.sx,floor(d1/2),floor(d2/2)]);
if isfield(options,'name') && ~isempty(options.name);
    name = [options.name,'_components'];
else
    name = [defoptions.name,'_components'];
end

% Count ROIs
nr = size(Ao,2);     % number of ROIs

% Get ROI masks
roi_masks=cell(1,length(nr));
for i=1:nr
    A_temp = full(reshape(Ao(:,i),d1,d2));
    A_temp = medfilt2(A_temp,[3,3]);
    A_temp = A_temp(:);
    roi_masks{i}=reshape(A_temp,d1,d2);
end

end