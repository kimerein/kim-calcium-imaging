clear;

curd=pwd;

cd('cvx');
cvx_setup;

cd(curd);
addpath(genpath('ca_source_extraction'));
addpath(genpath('utilities'));

% load file          
nam = 'combined1.tif';          % insert path to tiff stack here
sframe=1;						% user input: first frame to read (optional, default 1)
num2read=100000;					% user input: how many frames to read   (optional, default until the end)

Y = bigread2(nam,sframe,num2read);
Y = Y - min(Y(:)); 
if ~isa(Y,'double');    Y = double(Y);  end         % convert to double

[d1,d2,T] = size(Y);                                % dimensions of dataset
d = d1*d2;                                          % total number of pixels

% Set parameters

K = 200;                                           % number of components to be found
tau = 5;                                          % std of gaussian kernel (size of neuron) 
p = 1;                                            % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.75;                                  % merging threshold


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
% Data pre-processing

[P,Y] = preprocess_data(Y,p);

Cn =  reshape(P.sn,d1,d2);

% fast initialization of spatial components using greedyROI and HALS

[Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options,P);  % initialize

% update spatial components
Yr = reshape(Y,d,T);
clear Y;
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,Ain,P,options);

% update temporal components
P.p = 0;    % set AR temporarily to zero for speed
[C,f,P,S] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

% merge found components
[Am,Cm,K_m,merged_ROIs,P,Sm] = merge_components(Yr,A,b,C,f,P,S,options);

% repeat
P.p = p;    % restore AR value
[A2,b2,Cm] = update_spatial_components(Yr,Cm,f,Am,P,options);
[C2,f2,P,S2] = update_temporal_components(Yr,A2,b2,Cm,f,P,options);

% extract some useful values
[Ao,C_or,S_or,P] = order_ROIs(A2,C2,S2,P); % order components
K_m = size(C_or,1);
[Cf,Df] = extract_DF_F(Yr,[Ao,b2],[C_or;f2],K_m+1); % extract DF/F values 

nA = sqrt(sum(Ao.^2))';
[K,~] = size(C_or);
A_components = Ao/spdiags(nA,0,K,K);    % normalize spatial components to unit energy
C_components = spdiags(nA,0,K,K)*C_or;
Yk = (A_components'*Yr- (A_components'*A_components)*C_components - (A_components'*full(b2))*f2) + C_components;

% Raw trace for each component i is then Yk(i,:)/Df(i)
% Inferred trace for each component i is then C_components(i,:)/Df(i) or,
% equivalently, Cf

% output

% fprintf(1,'%6s\n','Yr');
% for i=1:size(Yr,1)
%     fprintf(1,'% 10.1f',Yr(i,:));
%     fprintf(1,'\n');
% end
% fprintf(1,'\n');

b2=b2';
fprintf(1,'%6s\n','b2');
for i=1:size(b2,1)
    fprintf(1,'% 10.2f',b2(i,:));
    fprintf(1,'\n');
end
fprintf(1,'\n');

fprintf(1,'%6s\n','f2');
for i=1:size(f2,1)
    fprintf(1,'% 6.6f',f2(i,:));
    fprintf(1,'\n');
end
fprintf(1,'\n');

fprintf(1,'%6s\n','Cn');
for i=1:size(Cn,1)
    fprintf(1,'% 10.4f',Cn(i,:));
    fprintf(1,'\n');
end
fprintf(1,'\n');

fprintf(1,'%6s\n','Yk');
for i=1:size(Yk,1)
    fprintf(1,'% 10.2f',Yk(i,:));
    fprintf(1,'\n');
end
fprintf(1,'\n');

fprintf(1,'%6s\n','Cf');
for i=1:size(Cf,1)
    fprintf(1,'% 10.4f',Cf(i,:));
    fprintf(1,'\n');
end
fprintf(1,'\n');

fprintf(1,'%6s\n','Df');
for i=1:size(Df,1)
    fprintf(1,'% 10.2f',Df(i,:));
    fprintf(1,'\n');
end
fprintf(1,'\n');

Ao=full(Ao');
fprintf(1,'%6s\n','Ao');
for i=1:size(Ao,1)
    fprintf(1,'% 10.4f',Ao(i,:));
    fprintf(1,'\n');
end
fprintf(1,'\n');

exit