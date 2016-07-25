function plot_dFoverF(C_df,acq_obj,nameOptoCommand,Yr,A_or,C_or,b2,f2,Cn,options,nameBehaviorCommand,saveDir,readin_data)
% C_df is deltaF over F in the form of a matrix where rows are neurons and
% columns are time points
% Last row of C_df is background
% acq_obj is Acquisition2P object from Harvey lab motion correction
% associated with this data set
% nameOptoCommand is string giving name associated with files containing
% opto stim data (e.g., 'Opto_Stim')
% Opto data assumed to be in the same folder as the movies originally used
% to create acq_obj, which is the same folder that contains the shutter
% data
% Yr,A_or,C_or,b2,f2,Cn,options are arguments to plot_components_GUI

if ~isempty(readin_data)
    beh_shutterTimesRemoved=readin_data.beh_shutterTimesRemoved;
    optoStimTypes=readin_data.optoStimTypes;
    optoMapping=readin_data.optoMapping;
    opto_shutterTimesRemoved=readin_data.opto_shutterTimesRemoved;
    if isfield(readin_data,'useComponents')
        useComponents=readin_data.useComponents;
    else
        useComponents=[];
    end
else
    % Read in opto stim data
    % If shutter data associated with acq_obj has already been saved, load it
    movName=acq_obj.Movies{1};
    dirBreaks=regexp(movName,'\','start');
    shutterPath=movName(1:dirBreaks(end));
    if ~isfield(acq_obj.sabaMetadata,'saveShutterDataFolder')
        shutterData=[];
    else
        listing=dir([shutterPath acq_obj.sabaMetadata.saveShutterDataFolder '\shutterData.mat']);
        if ~isempty(listing)
            a=load([shutterPath acq_obj.sabaMetadata.saveShutterDataFolder '\shutterData.mat']);
            shutterData=a.shutterData;
        else
            shutterData=[];
        end
    end
    [opto_shutterTimesRemoved,optoMapping,optoStimTypes]=loadOptoData(acq_obj,nameOptoCommand,shutterData);
    samplingRate=acq_obj.sabaMetadata.phys.settings.inputRate; % Get sampling rate of opto data
    times=0:1/samplingRate:(1/samplingRate)*length(optoMapping{1,1})-(1/samplingRate);
    [optoMapping,optoStimTypes]=groupOptoStimsGUI(optoMapping,optoStimTypes,times); % User decides how to group different opto stims for analysis
    beh_shutterTimesRemoved=loadGenericPhysData(acq_obj,nameBehaviorCommand,shutterData);
    % Time points when movie was shuttered have now been removed from opto data
    save([saveDir '\beh_shutterTimesRemoved.mat'],'beh_shutterTimesRemoved');
    save([saveDir '\optoStimTypes.mat'],'optoStimTypes');
    save([saveDir '\optoMapping.mat'],'optoMapping');
    save([saveDir '\opto_shutterTimesRemoved.mat'],'opto_shutterTimesRemoved');
end

% Now compare opto stim to dFoverF
% Use saved data (from motion correction) about which frames were removed
% from movies (when shutter on) to align opto, behavior and Ca2+ traces
[opto_stim,optoMovieLength]=alignToMovieFrames(opto_shutterTimesRemoved,C_df,acq_obj);
for i=1:length(beh_shutterTimesRemoved)
    vel=angular_velocity(beh_shutterTimesRemoved{i});
    vel=vel';
    vel=[vel*1000 vel(end)*1000]; % Pad at end so length of velocity vector matches length of position vector
    beh_shutterTimesRemoved{i}=vel;
end
[beh,behMovieLength]=alignToMovieFrames(beh_shutterTimesRemoved,C_df,acq_obj);

frameDuration=(acq_obj.sabaMetadata.acq.msPerLine/1000)*acq_obj.sabaMetadata.acq.linesPerFrame;
times=0:frameDuration:(size(C_df,2)-1)*frameDuration;

% Plot dFoverF GUI with opto data
plot_components_GUI_withopto_andBeh(Yr,A_or,C_or,b2,f2,Cn,options,opto_stim,beh,times);

% Get average opto-triggered responses
[optoTriggeredResponses,avOpto,trialByTrialBeh]=getOptoResponseAcrossCells(optoMovieLength,C_df,behMovieLength);

% Save progress to saveDir
savePartialProgress=1; % if savePartialProgress equals 1, save variables at this stage to saveDir
if savePartialProgress==1 && ~isempty(saveDir)
    if ~exist([saveDir '\partwayData_moviematched'],'dir')
        mkdir([saveDir '\partwayData_moviematched']);
    end
    save([saveDir '\partwayData_moviematched\trialByTrialBeh.mat'],'trialByTrialBeh');
    save([saveDir '\partwayData_moviematched\avOpto.mat'],'avOpto');
    save([saveDir '\partwayData_moviematched\optoMapping.mat'],'optoMapping');
    save([saveDir '\partwayData_moviematched\optoStimTypes.mat'],'optoStimTypes');
    save([saveDir '\partwayData_moviematched\acq_obj.mat'],'acq_obj');
    save([saveDir '\partwayData_moviematched\optoTriggeredResponses.mat'],'optoTriggeredResponses');
    if isfield(readin_data,'useComponents')
        save([saveDir '\partwayData_moviematched\useComponents.mat'],'useComponents');
    end 
end

% Find trials with desired behavior profile
profile=chooseBehaviorProfile(trialByTrialBeh,avOpto);

% Find trials with desired opto stim type
optoProfile=chooseOptoProfile(optoMapping,optoStimTypes,acq_obj);

% Combine behavioral and opto trial types
profile=(profile==1) & (optoProfile==1);

% Plot GUI with average opto-triggered responses from trials matching
% behavior profile
avResponses=takeTrialsForEachCell(optoTriggeredResponses,profile);
plot_avOptoTriggered_components_GUI(Yr,A_or,C_or,b2,f2,Cn,options,avResponses,avOpto,nanmean(trialByTrialBeh(profile,:),1),times);

% Plot average response across all cells
% Let user choose which cells to include
if ~isfield(readin_data,'useComponents')
    useComponents=chooseComponents(C_df);
end
% Get average
cellByCellAv=zeros(length(useComponents),size(avResponses{1},2));
for i=1:length(useComponents)
    cellByCellAv(i,:)=nanmean(avResponses{useComponents(i)},1);
end
figure();
hax=axes();
plotWStderr(hax,times,cellByCellAv,'k');
title('Average and Std Err across Cells');    

end

function outResponses=takeTrialsForEachCell(responses,useTrials)

outResponses=cell(1,length(responses));
for i=1:length(responses)
    r=responses{i};
    outR=r(useTrials,:);
    outResponses{i}=outR;
end

end

function [optoTriggeredResponses,avOpto,trialByTrialBeh]=getOptoResponseAcrossCells(optoMovieLength,C_df,behMovieLength)

optoTriggeredResponses=cell(1,size(C_df,1));
allOpto=nan(size(C_df,1),length(optoMovieLength{1}+10));
for i=1:size(C_df,1)
    % Consider each neuron's trace
    [optoTriggeredResponses{i},acrossTrialsOpto]=optoTriggeredResponse(optoMovieLength,C_df(i,:));
    avOpto=nanmean(acrossTrialsOpto,1);
    allOpto(i,1:length(avOpto))=avOpto;
end
beh=[];
for i=1:length(behMovieLength)
    beh=[beh behMovieLength{i}];
end
trialByTrialBeh=optoTriggeredResponse(optoMovieLength,beh);
avOpto=nanmean(allOpto,1);
avOpto=avOpto(~isnan(avOpto));

end

function [acrossTrialsResponse,acrossTrialsOpto]=optoTriggeredResponse(opto_stim,C_df)
% C_df is deltaFoverF trace from one neuron (vector)

acrossTrialsResponse=zeros(length(opto_stim),length(opto_stim{1})+10);
acrossTrialsOpto=zeros(length(opto_stim),length(opto_stim{1})+10);
j=1;
for i=1:length(opto_stim)
    currOpto=opto_stim{i};
    currCa=C_df(j:j+length(currOpto)-1);
    j=j+length(currOpto);
    if i==1
        optoAt=find(currOpto>0,1,'first');
        baseInds=optoAt-1;
    end
    currOptoAt=find(currOpto>0,1,'first');
    currBaseInds=currOptoAt-1;
    if currBaseInds<baseInds
        acrossTrialsResponse(i,:)=[zeros(1,baseInds-currBaseInds) currCa(1:end-(baseInds-currBaseInds))];
        acrossTrialsOpto(i,:)=[zeros(1,baseInds-currBaseInds) currOpto(1:end-(baseInds-currBaseInds))];
    elseif currBaseInds>baseInds
        acrossTrialsResponse(i,:)=[currCa(currBaseInds-baseInds+1:end) zeros(1,currBaseInds-baseInds)];
        acrossTrialsOpto(i,:)=[currOpto(currBaseInds-baseInds+1:end) zeros(1,currBaseInds-baseInds)];
    else
        acrossTrialsResponse(i,:)=currCa;
        acrossTrialsOpto(i,:)=currOpto;
    end
end
    
end

function newSignal=putDeltasIntoResample(signal,newSize,deltaThresh)

times_signal=linspace(0,100,length(signal));
times_newSize=linspace(0,100,newSize);

% Find closest times in times_newSize when signal exceeds deltaThresh
isHigh=zeros(1,newSize);
timesHigh=times_signal(signal>deltaThresh);
valHigh=signal(signal>deltaThresh);
for i=1:length(timesHigh)
    temp=abs(times_newSize-timesHigh(i));
    [~,idx]=min(temp);
%     isHigh(idx)=1;
    isHigh(idx)=valHigh(i);
end

newSignal=zeros(1,newSize);
% newSignal(isHigh==1)=1;
newSignal(isHigh>0)=isHigh(isHigh>0);

end

function plot_avOptoTriggered_components_GUI(Y,A,C,b,f,Cn,options,acrossTrials,opto,beh,times)
defoptions = CNMFSetParms;
if nargin < 7 || isempty(options); options = []; end
if ~isfield(options,'d1') || isempty(options.d1); d1 = input('What is the total number of rows? \n'); else d1 = options.d1; end          % # of rows
if ~isfield(options,'d2') || isempty(options.d2); d2 = input('What is the total number of columns? \n'); else d2 = options.d2; end          % # of columns
% if ~isfield(options,'normalize') || isempty(options.normalize); options.normalize = ones(size(A,1),1); end
%     sn = options.normalize;
if ~isfield(options,'plot_df') || isempty(options.plot_df); options.df = defoptions.plot_df; end
plot_df = options.plot_df;
if ~isfield(options,'make_gif') || isempty(options.make_gif); options.make_gif = defoptions.make_gif; end
make_gif = options.make_gif;
if ~isfield(options,'save_avi') || isempty(options.save_avi); options.save_avi = defoptions.save_avi; end
save_avi = options.save_avi;
if ~isfield(options,'sx') || isempty(options.sx); options.sx = defoptions.sx; end
sx = min([options.sx,floor(d1/2),floor(d2/2)]);
%if ~isfield(options,'pause_time') || isempty(options.pause_time); options.pause_time = defoptions.pause_time; end
%pause_time = options.pause_time;
if isfield(options,'name') && ~isempty(options.name);
    name = [options.name,'_components'];
else
    name = [defoptions.name,'_components'];
end

T = size(C,2);
if ndims(Y) == 3
    Y = reshape(Y,d1*d2,T);
end
if nargin < 6 || isempty(Cn);
    Cn = reshape(mean(Y,2),d1,d2);
end

nA = sqrt(sum(A.^2))';
[K,~] = size(C);
A = A/spdiags(nA,0,K,K);    % normalize spatial components to unit energy
C = spdiags(nA,0,K,K)*C;

nr = size(A,2);     % number of ROIs
nb = size(f,1);     % number of background components
%nA = full(sum(A.^2))';  % energy of each row
%Y_r = spdiags(nA,0,nr,nr)\(A'*Y- (A'*A)*C - (A'*full(b))*f) + C; 
Y_r = (A'*Y- (A'*A)*C - (A'*full(b))*f) + C;

if plot_df
    [~,Df] = extract_DF_F(Y,[A,b],[C;f],size(A,2)+1);
else
    Df = ones(size(A,2)+1,1);
end

if save_avi
    vidObj = VideoWriter([name,'.avi']);
    set(vidObj,'FrameRate',1);
    open(vidObj);
end
thr = 0.95;
fig = figure('Visible','off');
set(gcf,'Position',2*[300,300,960,480]);
set(gcf,'PaperPosition',2*[300,300,960,480]);
int_x = zeros(nr,2*sx);
int_y = zeros(nr,2*sx);
cm = com(A,d1,d2);

% Create a figure and axes
% ax = axes('Units','DF/F');

% Create slider
sld = uicontrol('Style', 'slider',...
    'Min',1,'Max',nr+nb,'Value',1,'SliderStep',[1/(nr+nb-1) 1],...
    'Position', [150 20 800 20],...
    'Callback', @surfzlim);

% Add a text uicontrol to label the slider.
txt = uicontrol('Style','text',...
    'Position',[400 45 120 20],...
    'String','Component');

% Make figure visble after adding all components
fig.Visible = 'on';
plot_component(1)

% This code uses dot notation to set properties.
% Dot notation runs in R2014b and later.
% For R2014a and earlier: set(f,'Visible','on');


    function surfzlim(source,callbackdata)
        i = source.Value;
        plot_component(round(i))
        % For R2014a and earlier:
        % i = get(source,'Value');
        
        
%         if save_avi
%             currFrame = getframe(fig);
%             writeVideo(vidObj,currFrame);
%         else
%             pause(0.05);
%         end
    end

    function plot_component(i)
       if i <= nr
            subplot(3,2,5);
            Atemp = reshape(A(:,i),d1,d2);
            int_x(i,:) = round(cm(i,1)) + (-(sx-1):sx);
            if int_x(i,1)<1
                int_x(i,:) = int_x(i,:) + 1 - int_x(i,1);
            end
            if int_x(i,end)>d1
                int_x(i,:) = int_x(i,:) - (int_x(i,end)-d1);
            end
            int_y(i,:) = round(cm(i,2)) + (-(sx-1):sx);
            if int_y(i,1)<1
                int_y(i,:) = int_y(i,:) + 1 - int_y(i,1);
            end
            if int_y(i,end)>d2
                int_y(i,:) = int_y(i,:) - (int_y(i,end)-d2);
            end
            Atemp = Atemp(int_x(i,:),int_y(i,:));
            imagesc(int_x(i,:),int_y(i,:),Atemp); axis square;
        end
        subplot(3,2,[1,3]);
        if i <= nr
            cla
            imagesc(2*Cn); axis equal; axis tight; axis off; hold on;
            A_temp = full(reshape(A(:,i),d1,d2));
            A_temp = medfilt2(A_temp,[3,3]);
            A_temp = A_temp(:);
            [temp,ind] = sort(A_temp(:).^2,'ascend');
            temp =  cumsum(temp);
            ff = find(temp > (1-thr)*temp(end),1,'first');
            if ~isempty(ff)
                [~,ww] = contour(reshape(A_temp,d1,d2),[0,0]+A_temp(ind(ff)),'LineColor','k');
                ww.LineWidth = 2;
            end
            title(sprintf('Component %i ',i),'fontsize',16,'fontweight','bold'); drawnow; %pause;
        else
            cla
            imagesc(reshape(b(:,i-nr),d1,d2)); axis equal; axis tight;
            title('Background component','fontsize',16,'fontweight','bold'); drawnow;
        end
        subplot(3,2,[2,4,6]);
        if i <= nr
%             plot(1:T,Y_r(i,:)/Df(i),'linewidth',2); hold all; 
            currAv=acrossTrials{i};
            plot(times,nanmean(currAv,1),'linewidth',2,'Color','r');
            hold all;
%             Overlay opto stim times
            plot(times,(opto./max(opto))*max(nanmean(currAv,1)),'linewidth',2);
            plot(times,beh,'linewidth',1);
            if plot_df
                title(sprintf('Component %i (calcium DF/F value)',i),'fontsize',16,'fontweight','bold');
            else
                title(sprintf('Component %i (calcium raw value)',i),'fontsize',16,'fontweight','bold');
            end
            leg = legend('Inferred','Opto Stim','Angular Velocity');
            set(leg,'FontSize',14,'FontWeight','bold'); 
            plot(times,nanmean(currAv,1)+nanstd(currAv,[],1)./sqrt(size(currAv,1)),'linewidth',1,'Color','r');
            plot(times,nanmean(currAv,1)-nanstd(currAv,[],1)./sqrt(size(currAv,1)),'linewidth',1,'Color','r');
            xlabel('Time (s)');
            drawnow;
            hold off;
            if make_gif
%                 frame = getframe(fig); %getframe(1);
%                 im = frame2im(frame);
%                 [imind,clm] = rgb2ind(im,256);
%                 if i == 1;
%                     imwrite(imind,clm,[name,'.gif'],'gif', 'Loopcount',inf);
%                 else
%                     imwrite(imind,clm,[name,'.gif'],'gif','WriteMode','append');
%                 end
            else
%                 if i < nr+nb && ~save_avi
%                     fprintf('component %i. Press any key to continue.. \n', i);
%                     if pause_time == Inf;
%                         pause;
%                     else
%                         pause(pause_time);
%                     end

%                 end
            end
        else
            plot(1:T,f(i-nr,:)); title('Background activity','fontsize',16,'fontweight','bold');
            drawnow;
            if make_gif
                frame = getframe(fig); %getframe(1);
                im = frame2im(frame);
                [imind,clm] = rgb2ind(im,256);
                if i == 1;
                    imwrite(imind,clm,[name,'.gif'],'gif', 'Loopcount',inf);
                else
                    imwrite(imind,clm,[name,'.gif'],'gif','WriteMode','append');
                end
            else
%                 if i < nr+nb && ~save_avi
%                     fprintf('background component %i. Press any key to continue.. \n', i-nr);
% %                     if pause_time == Inf;
% %                         pause;
% %                     else
% %                         pause(pause_time);%                     end
%                 end
            end
        end 
    end
end

function plot_components_GUI_withopto_andBeh(Y,A,C,b,f,Cn,options,opto_stim,beh_resample,times)
defoptions = CNMFSetParms;
if nargin < 7 || isempty(options); options = []; end
if ~isfield(options,'d1') || isempty(options.d1); d1 = input('What is the total number of rows? \n'); else d1 = options.d1; end          % # of rows
if ~isfield(options,'d2') || isempty(options.d2); d2 = input('What is the total number of columns? \n'); else d2 = options.d2; end          % # of columns
% if ~isfield(options,'normalize') || isempty(options.normalize); options.normalize = ones(size(A,1),1); end
%     sn = options.normalize;
if ~isfield(options,'plot_df') || isempty(options.plot_df); options.df = defoptions.plot_df; end
plot_df = options.plot_df;
if ~isfield(options,'make_gif') || isempty(options.make_gif); options.make_gif = defoptions.make_gif; end
make_gif = options.make_gif;
if ~isfield(options,'save_avi') || isempty(options.save_avi); options.save_avi = defoptions.save_avi; end
save_avi = options.save_avi;
if ~isfield(options,'sx') || isempty(options.sx); options.sx = defoptions.sx; end
sx = min([options.sx,floor(d1/2),floor(d2/2)]);
%if ~isfield(options,'pause_time') || isempty(options.pause_time); options.pause_time = defoptions.pause_time; end
%pause_time = options.pause_time;
if isfield(options,'name') && ~isempty(options.name);
    name = [options.name,'_components'];
else
    name = [defoptions.name,'_components'];
end

T = size(C,2);
if ndims(Y) == 3
    Y = reshape(Y,d1*d2,T);
end
if nargin < 6 || isempty(Cn);
    Cn = reshape(mean(Y,2),d1,d2);
end

nA = sqrt(sum(A.^2))';
[K,~] = size(C);
A = A/spdiags(nA,0,K,K);    % normalize spatial components to unit energy
C = spdiags(nA,0,K,K)*C;

nr = size(A,2);     % number of ROIs
nb = size(f,1);     % number of background components
%nA = full(sum(A.^2))';  % energy of each row
%Y_r = spdiags(nA,0,nr,nr)\(A'*Y- (A'*A)*C - (A'*full(b))*f) + C; 
Y_r = (A'*Y- (A'*A)*C - (A'*full(b))*f) + C;

if plot_df
    [~,Df] = extract_DF_F(Y,[A,b],[C;f],size(A,2)+1);
else
    Df = ones(size(A,2)+1,1);
end

if save_avi
    vidObj = VideoWriter([name,'.avi']);
    set(vidObj,'FrameRate',1);
    open(vidObj);
end
thr = 0.95;
fig = figure('Visible','off');
set(gcf,'Position',2*[300,300,960,480]);
set(gcf,'PaperPosition',2*[300,300,960,480]);
int_x = zeros(nr,2*sx);
int_y = zeros(nr,2*sx);
cm = com(A,d1,d2);

% Create a figure and axes
% ax = axes('Units','DF/F');

% Create slider
sld = uicontrol('Style', 'slider',...
    'Min',1,'Max',nr+nb,'Value',1,'SliderStep',[1/(nr+nb-1) 1],...
    'Position', [150 20 800 20],...
    'Callback', @surfzlim);

% Add a text uicontrol to label the slider.
txt = uicontrol('Style','text',...
    'Position',[400 45 120 20],...
    'String','Component');

% Make figure visble after adding all components
fig.Visible = 'on';
plot_component(1)

% This code uses dot notation to set properties.
% Dot notation runs in R2014b and later.
% For R2014a and earlier: set(f,'Visible','on');


    function surfzlim(source,callbackdata)
        i = source.Value;
        plot_component(round(i))
        % For R2014a and earlier:
        % i = get(source,'Value');
        
        
%         if save_avi
%             currFrame = getframe(fig);
%             writeVideo(vidObj,currFrame);
%         else
%             pause(0.05);
%         end
    end

    function plot_component(i)
       if i <= nr
            subplot(3,2,5);
            Atemp = reshape(A(:,i),d1,d2);
            int_x(i,:) = round(cm(i,1)) + (-(sx-1):sx);
            if int_x(i,1)<1
                int_x(i,:) = int_x(i,:) + 1 - int_x(i,1);
            end
            if int_x(i,end)>d1
                int_x(i,:) = int_x(i,:) - (int_x(i,end)-d1);
            end
            int_y(i,:) = round(cm(i,2)) + (-(sx-1):sx);
            if int_y(i,1)<1
                int_y(i,:) = int_y(i,:) + 1 - int_y(i,1);
            end
            if int_y(i,end)>d2
                int_y(i,:) = int_y(i,:) - (int_y(i,end)-d2);
            end
            Atemp = Atemp(int_x(i,:),int_y(i,:));
            imagesc(int_x(i,:),int_y(i,:),Atemp); axis square;
        end
        subplot(3,2,[1,3]);
        if i <= nr
            cla
            imagesc(2*Cn); axis equal; axis tight; axis off; hold on;
            A_temp = full(reshape(A(:,i),d1,d2));
            A_temp = medfilt2(A_temp,[3,3]);
            A_temp = A_temp(:);
            [temp,ind] = sort(A_temp(:).^2,'ascend');
            temp =  cumsum(temp);
            ff = find(temp > (1-thr)*temp(end),1,'first');
            if ~isempty(ff)
                [~,ww] = contour(reshape(A_temp,d1,d2),[0,0]+A_temp(ind(ff)),'LineColor','k');
                ww.LineWidth = 2;
            end
            title(sprintf('Component %i ',i),'fontsize',16,'fontweight','bold'); drawnow; %pause;
        else
            cla
            imagesc(reshape(b(:,i-nr),d1,d2)); axis equal; axis tight;
            title('Background component','fontsize',16,'fontweight','bold'); drawnow;
        end
        subplot(3,2,[2,4,6]);
        if i <= nr
            plot(times,Y_r(i,:)/Df(i),'linewidth',2); 
            hold all;
            plot(times,C(i,:)/Df(i),'linewidth',2);
            plot(times,beh_resample,'linewidth',1);
            plot(times,opto_stim,'linewidth',2);
            xlabel('Time (s)');
            % Overlay opto stim times
            if plot_df
                title(sprintf('Component %i (calcium DF/F value)',i),'fontsize',16,'fontweight','bold');
            else
                title(sprintf('Component %i (calcium raw value)',i),'fontsize',16,'fontweight','bold');
            end
            leg = legend('Raw trace (filtered)','Inferred','Angular Velocity','Opto Stim');
            set(leg,'FontSize',14,'FontWeight','bold');
            drawnow;
            hold off;
            if make_gif
%                 frame = getframe(fig); %getframe(1);
%                 im = frame2im(frame);
%                 [imind,clm] = rgb2ind(im,256);
%                 if i == 1;
%                     imwrite(imind,clm,[name,'.gif'],'gif', 'Loopcount',inf);
%                 else
%                     imwrite(imind,clm,[name,'.gif'],'gif','WriteMode','append');
%                 end
            else
%                 if i < nr+nb && ~save_avi
%                     fprintf('component %i. Press any key to continue.. \n', i);
%                     if pause_time == Inf;
%                         pause;
%                     else
%                         pause(pause_time);
%                     end

%                 end
            end
        else
            plot(1:T,f(i-nr,:)); title('Background activity','fontsize',16,'fontweight','bold');
            drawnow;
            if make_gif
                frame = getframe(fig); %getframe(1);
                im = frame2im(frame);
                [imind,clm] = rgb2ind(im,256);
                if i == 1;
                    imwrite(imind,clm,[name,'.gif'],'gif', 'Loopcount',inf);
                else
                    imwrite(imind,clm,[name,'.gif'],'gif','WriteMode','append');
                end
            else
%                 if i < nr+nb && ~save_avi
%                     fprintf('background component %i. Press any key to continue.. \n', i-nr);
% %                     if pause_time == Inf;
% %                         pause;
% %                     else
% %                         pause(pause_time);%                     end
%                 end
            end
        end 
    end
end

