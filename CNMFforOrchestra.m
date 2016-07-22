clear;
%% load file

addpath(genpath('utilities'));
             
nam = 'combined1.tif';          % insert path to tiff stack here
sframe=1;						% user input: first frame to read (optional, default 1)
num2read=100000;					% user input: how many frames to read   (optional, default until the end)

Y = bigread2(nam,sframe,num2read);
Y = Y - min(Y(:)); 
if ~isa(Y,'double');    Y = double(Y);  end         % convert to double

[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;                                          % total number of pixels

%% Set parameters

K = 50;                                           % number of components to be found
tau = 5;                                          % std of gaussian kernel (size of neuron) 
p = 1;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.75;                                  % merging threshold
saveDir='outputOfCNMF';

if ~exists(saveDir)
    status=mkdir(saveDir);
end

options = CNMFSetParms(...                      
    'd1',d1,'d2',d2,...                         % dimensions of datasets
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
%% Data pre-processing

[P,Y] = preprocess_data(Y,p);

%% fast initialization of spatial components using greedyROI and HALS

[Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options,P);  % initialize

%% update spatial components
Yr = reshape(Y,d,T);
clear Y;
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,Ain,P,options);

%% update temporal components
P.p = 0;    % set AR temporarily to zero for speed
[C,f,P,S] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

%% merge found components
[Am,Cm,K_m,merged_ROIs,P,Sm] = merge_components(Yr,A,b,C,f,P,S,options);

%% repeat
P.p = p;    % restore AR value
[A2,b2,Cm] = update_spatial_components(Yr,Cm,f,Am,P,options);
[C2,f2,P,S2] = update_temporal_components(Yr,A2,b2,Cm,f,P,options);

%% output
[A_or,C_or,S_or,P] = order_ROIs(A2,C2,S2,P); % order components
K_m = size(C_or,1);
[C_df,~] = extract_DF_F(Yr,[A_or,b2],[C_or;f2],K_m+1); % extract DF/F values (optional)