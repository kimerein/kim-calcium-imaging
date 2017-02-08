function dFoverF_viewer(C_df,acq_obj,nameOptoCommand,Yk,A_or,Cn,b2,f2,Df,options,nameBehaviorCommand,saveDir,readin_data)
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
% Yk,A_or,Cn,b2,f2,Df,options are arguments to plot_components_GUI

getOptoMapping=1;
getComponents=1;
if ~isempty(readin_data)
    if isfield(readin_data,'useComponents')
        useComponents=readin_data.useComponents;
        getComponents=0;
    else
        useComponents=[];
        getComponents=1;
    end
    if isfield(readin_data,'optoMapping')
        optoStimTypes=readin_data.optoStimTypes;
        optoMapping=readin_data.optoMapping;
        getOptoMapping=0;
    else
        optoStimTypes=[];
        optoMapping=[];
        getOptoMapping=1;
    end
end
getData=1;
if ~isempty(readin_data) 
    if ~isfield(readin_data,'opto_shutterTimesRemoved')
        getData=1;
    else
        getData=0;
    end
else
    getData=1;
end
if getData==0
    beh_shutterTimesRemoved=readin_data.beh_shutterTimesRemoved;
    optoStimTypes=readin_data.optoStimTypes;
    optoMapping=readin_data.optoMapping;
    opto_shutterTimesRemoved=readin_data.opto_shutterTimesRemoved;    
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
    
    if getOptoMapping==1
        [optoMapping,optoStimTypes]=groupOptoStimsGUI(optoMapping,optoStimTypes,times); % User decides how to group different opto stims for analysis
    end
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
[opto_stim,optoMovieLength]=alignToMovieFrames(opto_shutterTimesRemoved,C_df,acq_obj,1);
for i=1:length(beh_shutterTimesRemoved)
    % Often problematic transients at beginning of phys recording, so zero
    % out
    temp=beh_shutterTimesRemoved{i};
    temp(1:500)=temp(500);
    temp(end-500:end)=temp(end-500);
    vel=angular_velocity(temp,acq_obj);
    vel=vel';
    vel=[vel*1000 vel(end)*1000]; % Pad at end so length of velocity vector matches length of position vector
    beh_shutterTimesRemoved{i}=vel; 
end 
[beh,behMovieLength]=alignToMovieFrames(beh_shutterTimesRemoved,C_df,acq_obj,0);

frameDuration=(acq_obj.sabaMetadata.acq.msPerLine/1000)*acq_obj.sabaMetadata.acq.linesPerFrame;
times=0:frameDuration:(size(C_df,2)-1)*frameDuration;

% Plot dFoverF GUI with opto data
[~,~,~,~,~,~,responsesType]=analysisSettings();
if strcmp(responsesType.type,'spikes')
    % Transform background-subtracted to spikes using constrained_foopsi deconvolution
    rawResponses=nan(size(C_df,1),size(Yk,2));
    for i=1:size(C_df,1)
        rawResponses(i,:)=Yk(i,:)/Df(i);
    end
    % Transform inferred to spiking using constrained_foopsi deconvolution
    temp=nan(size(C_df));
    for i=1:size(C_df,1)
        temp(i,:)=get_spikes_from_Ca2(C_df(i,:),options);
    end
    C_df=temp;
end
plot_components_GUI_withopto_andBeh(C_df,Cn,f2,A_or,b2,Yk,Df,options,opto_stim,beh,times);

% Get average opto-triggered responses 
if strcmp(responsesType.type,'raw')
    rawResponses=nan(size(C_df,1),size(Yk,2));
    for i=1:size(C_df,1)-1
        rawResponses(i,:)=Yk(i,:)/Df(i);
    end
    [optoTriggeredResponses,acrossTrialsOpto,trialByTrialBeh]=getOptoResponseAcrossCells(optoMovieLength,rawResponses,behMovieLength);
elseif strcmp(responsesType.type,'inferred')
    [optoTriggeredResponses,acrossTrialsOpto,trialByTrialBeh]=getOptoResponseAcrossCells(optoMovieLength,C_df,behMovieLength);
elseif strcmp(responsesType.type,'spikes')
%     [optoTriggeredResponses,acrossTrialsOpto,trialByTrialBeh]=getOptoResponseAcrossCells(optoMovieLength,C_df,behMovieLength);
    rawResponses=nan(size(C_df,1),size(Yk,2));
    for i=1:size(C_df,1)
        rawResponses(i,:)=Yk(i,:)/Df(i);
    end    
    [optoTriggeredResponses,acrossTrialsOpto,trialByTrialBeh]=getOptoResponseAcrossCells(optoMovieLength,rawResponses,behMovieLength);
else
    error('Do not recognize response.type in analysis settings');
end

if strcmp(responsesType.type,'spikes')
    for i=1:length(optoTriggeredResponses)
        curr=optoTriggeredResponses{i};
        stitched=[];
        for j=1:size(curr,1)
            if j==1
                stitched=[stitched curr(j,:)];
            else
                stitched=[stitched curr(j,:)-curr(j,1)+stitched(end)];
            end
        end
        optoTriggeredResponses{i}=get_spikes_from_Ca2(stitched,options);
    end
end
            

% Save progress to saveDir
savePartialProgress=1; % if savePartialProgress equals 1, save variables at this stage to saveDir
if savePartialProgress==1 && ~isempty(saveDir)
    if ~exist([saveDir '\partwayData_moviematched'],'dir')
        mkdir([saveDir '\partwayData_moviematched']);
    end
    save([saveDir '\partwayData_moviematched\trialByTrialBeh.mat'],'trialByTrialBeh');
    save([saveDir '\partwayData_moviematched\acrossTrialsOpto.mat'],'acrossTrialsOpto');
    save([saveDir '\partwayData_moviematched\optoMapping.mat'],'optoMapping');
    save([saveDir '\partwayData_moviematched\optoStimTypes.mat'],'optoStimTypes');
    save([saveDir '\partwayData_moviematched\acq_obj.mat'],'acq_obj');
    save([saveDir '\partwayData_moviematched\optoTriggeredResponses.mat'],'optoTriggeredResponses');
    if isfield(readin_data,'useComponents')
        save([saveDir '\partwayData_moviematched\useComponents.mat'],'useComponents');
    end 
end

% Find trials with desired behavior profile
profile=chooseBehaviorProfile(trialByTrialBeh,nanmean(acrossTrialsOpto,1));

% Find trials with desired opto stim type
optoProfile=chooseOptoProfile(optoMapping,optoStimTypes,acq_obj);

% Combine behavioral and opto trial types
profile=(profile==1) & (optoProfile==1);

% Plot GUI with average opto-triggered responses from trials matching
% behavior profile
avResponses=takeTrialsForEachCell(optoTriggeredResponses,profile);
times=0:frameDuration:(size(avResponses{1},2)-1)*frameDuration;
optoTriggered_viewer(C_df,Cn,f2,A_or,b2,options,avResponses,nanmean(acrossTrialsOpto(profile,:),1),nanmean(trialByTrialBeh(profile,:),1),times);

% Plot average response across all cells
% Let user choose which cells to include
if ~isfield(readin_data,'useComponents')
    useComponents=chooseComponents(C_df);
    save([saveDir '\partwayData_moviematched\useComponents.mat'],'useComponents');
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
hold on; 
plot(hax,times,nanmean(acrossTrialsOpto(profile,:),1),'Color','c');

end

function outResponses=takeTrialsForEachCell(responses,useTrials)

outResponses=cell(1,length(responses));
for i=1:length(responses)
    r=responses{i};
    outR=r(useTrials,:);
    outResponses{i}=outR;
end

end

function [optoTriggeredResponses,acrossTrialsOpto,trialByTrialBeh,optoAt]=getOptoResponseAcrossCells(optoMovieLength,C_df,behMovieLength)

baselineSubtract=0;

optoTriggeredResponses=cell(1,size(C_df,1));
acrossTrialsOpto=cell(1,size(C_df,1));
for i=1:size(C_df,1)
    % Consider each neuron's trace
    [temp,acrossTrialsOpto,baseInds]=optoTriggeredResponse(optoMovieLength,C_df(i,:));
    if baselineSubtract==1
        optoTriggeredResponses{i}=temp-repmat(nanmean(temp(:,baseInds),2),1,size(temp,2));
    else
        optoTriggeredResponses{i}=temp;
    end
end
beh=[];
for i=1:length(behMovieLength)
    beh=[beh behMovieLength{i}];
end
trialByTrialBeh=optoTriggeredResponse(optoMovieLength,beh);

end

function [acrossTrialsResponse,acrossTrialsOpto,baseInds]=optoTriggeredResponse(opto_stim,C_df)
% C_df is deltaFoverF trace from one neuron (vector)

[~,~,~,~,~,~,response]=analysisSettings;

j=1;
% acrossTrialsResponse=zeros(length(opto_stim),length(opto_stim{1})+10);
% acrossTrialsOpto=zeros(length(opto_stim),length(opto_stim{1})+10);
for i=1:length(opto_stim)
    if i==1
        prevOpto=nan(size(opto_stim{i}));
        prevCa=nan(size(C_df(j:j+length(prevOpto)-1)));
    else 
        prevOpto=currOpto;
        prevCa=currCa;
    end
    currOpto=opto_stim{i};
    currCa=C_df(j:j+length(currOpto)-1);
    j=j+length(currOpto);
    if i==1
        optoAt=find(currOpto>0,1,'first');
        baseInds=optoAt-1;
        if response.include_previous_trial==1
            acrossTrialsResponse=zeros(length(opto_stim),length(currCa)+response.n_previous_trial_inds);
            acrossTrialsOpto=zeros(length(opto_stim),length(currCa)+response.n_previous_trial_inds);
        else
            acrossTrialsResponse=zeros(length(opto_stim),length(currCa));
            acrossTrialsOpto=zeros(length(opto_stim),length(currCa));
        end
    end
    currOptoAt=find(currOpto>0,1,'first');
    currBaseInds=currOptoAt-1;
    if currBaseInds<baseInds
%         acrossTrialsResponse(i,:)=[zeros(1,baseInds-currBaseInds) currCa(1:end-(baseInds-currBaseInds))];
%         acrossTrialsOpto(i,:)=[zeros(1,baseInds-currBaseInds) currOpto(1:end-(baseInds-currBaseInds))];
        currCa=[zeros(1,baseInds-currBaseInds) currCa(1:end-(baseInds-currBaseInds))];
        currOpto=[zeros(1,baseInds-currBaseInds) currOpto(1:end-(baseInds-currBaseInds))];
    elseif currBaseInds>baseInds
%         acrossTrialsResponse(i,:)=[currCa(currBaseInds-baseInds+1:end) zeros(1,currBaseInds-baseInds)];
%         acrossTrialsOpto(i,:)=[currOpto(currBaseInds-baseInds+1:end) zeros(1,currBaseInds-baseInds)];
        currCa=[currCa(currBaseInds-baseInds+1:end) zeros(1,currBaseInds-baseInds)];
        currOpto=[currOpto(currBaseInds-baseInds+1:end) zeros(1,currBaseInds-baseInds)];
    end
    % Add some data from previous trial to beginning of this trial?
    if response.include_previous_trial==1
        currCa=[prevCa(end-response.n_previous_trial_inds:end) currCa];
        currOpto=[prevOpto(end-response.n_previous_trial_inds:end) currOpto];
    end
    if size(acrossTrialsResponse,2)>length(currCa) % This trial is shorter ... fill with nans
        temp=nan(1,size(acrossTrialsResponse,2));
        temp(1:length(currCa))=currCa;
        acrossTrialsResponse(i,:)=temp;
        temp=nan(1,size(acrossTrialsOpto,2));
        temp(1:length(currOpto))=currOpto;
        acrossTrialsOpto(i,:)=temp;
    elseif length(currCa)>size(acrossTrialsResponse,2) % This trial is longer ... truncate to match other trials' lengths
        acrossTrialsResponse(i,:)=currCa(1:size(acrossTrialsResponse,2));
        acrossTrialsOpto(i,:)=currOpto(1:size(acrossTrialsOpto,2));
    else % Lengths match
        acrossTrialsResponse(i,:)=currCa;
        acrossTrialsOpto(i,:)=currOpto;
    end
end
 
if response.include_previous_trial==1
    baseInds=baseInds+response.n_previous_trial_inds;
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
