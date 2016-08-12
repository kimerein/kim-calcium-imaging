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