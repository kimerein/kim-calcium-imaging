function [withinCellResponses,withinCellStats,withinCellAverages,times,optoForProfile,out,behForProfile]=analysisSecondHalf(loadDir,saveData)

% Load in data from first half of analysis
listing=dir(loadDir);
isMat=zeros(1,length(listing));
for i=1:length(listing)
    if ~isempty(regexp(listing(i).name,'.mat'))
        isMat(i)=1;
    end
end
listing=listing(logical(isMat));

for i=1:length(listing)
    load([loadDir '\' listing(i).name]);
end

% Begin rest of analysis

% Choose behavioral profiles and opto stim types
[bset,oset,cset,sorting]=analysisSettings();
% Iterate analysis on each combo

behProfiles=bset.profiles(bset.show_profiles==1);
optProfiles=oset.profiles(oset.show_profiles==1);
    
% Times
frameDuration=(acq_obj.sabaMetadata.acq.msPerLine/1000)*acq_obj.sabaMetadata.acq.linesPerFrame;
times=0:frameDuration:(size(acrossTrialsOpto,2)-1)*frameDuration;

withinCellResponses=cell(length(behProfiles),length(optProfiles));
withinCellAverages=cell(length(behProfiles),length(optProfiles));
withinCellChanges=cell(length(behProfiles),length(optProfiles));
withinCellStats=cell(length(behProfiles),length(optProfiles));
withinCellPlotOutput=cell(length(behProfiles),length(optProfiles));
optoForProfile=cell(length(behProfiles),length(optProfiles));
behForProfile=cell(length(behProfiles),length(optProfiles));
usedProfiles=cell(1,length(behProfiles)*length(optProfiles));
countProfiles=1;
for i=1:length(behProfiles)
    % Find trials with desired behavior profile  
    behOnlyProfile=getBehaviorProfile(trialByTrialBeh,nanmean(acrossTrialsOpto,1)./max(nanmean(acrossTrialsOpto,1)),behProfiles,i);
    for j=1:length(optProfiles)
        % Find trials with desired opto stim type
        optoOnlyProfile=optoStimProfile(optoStimTypes,optProfiles{j});
        % Combine behavioral and opto trial types
        profile=(behOnlyProfile==1) & (optoOnlyProfile==1);
        usedProfiles{countProfiles}=profile;
        countProfiles=countProfiles+1;
        % Get response for each cell over trials matching profile
        responses=takeTrialsForEachCell(optoTriggeredResponses,profile);
        optoForProfile{i,j}=acrossTrialsOpto(profile,:);
        withinCellResponses{i,j}=responses(useComponents);
        
        % Get average behavior for this profile
        behForProfile{i,j}=nanmean(trialByTrialBeh(profile,:),1);
        
        % Get average response across all cells
        if sum(profile)==0
            % No trials match this profile
            cellByCellAv=zeros(length(useComponents),size(optoTriggeredResponses{1},2));
        else
            cellByCellAv=zeros(length(useComponents),size(responses{1},2));
            for k=1:length(useComponents)
                cellByCellAv(k,:)=nanmean(responses{useComponents(k)},1);
            end
        end
        withinCellAverages{i,j}=cellByCellAv;
    
        % Get change in response during timewindow
        if cset.baseline_subtract==1
            withinCellChanges{i,j}=nanmean(cellByCellAv(:,times>=cset.timewindow(1) & times<=cset.timewindow(2)),2)-nanmean(cellByCellAv(:,times>=cset.baselinewindow(1) & times<=cset.baselinewindow(2)),2);
        else
            withinCellChanges{i,j}=nanmean(cellByCellAv(:,times>=cset.timewindow(1) & times<=cset.timewindow(2)),2);
        end
        
        % Get statistics on change in response during timewindow
        cellByCellStat=zeros(length(useComponents),1);
        responses=withinCellResponses{i,j};
        for k=1:length(responses)
            currResponse=responses{k};
            currTimewindow=nanmean(currResponse(:,times>=cset.timewindow(1) & times<=cset.timewindow(2)),2);
            currBaselinewindow=nanmean(currResponse(:,times>=cset.baselinewindow(1) & times<=cset.baselinewindow(2)),2);
            switch cset.stats.test_type
                case 'signrank'
                    % Paired non-parametric
%                     if cset.stats.vs_baseline==1
%                         if isempty(currTimewindow-currBaselinewindow)
%                             p=nan;
%                         else
%                             p=signrank(currTimewindow-currBaselinewindow);
%                         end
%                     else
%                         currOtherwindow=nanmean(currResponse(:,times>=cset.stats.vs_othertimewindow(1) & times<=cset.stats.vs_othertimewindow(2)),2);
%                         if isempty(currTimewindow-currOtherwindow)
%                             p=nan;
%                         else
%                             p=signrank(currTimewindow-currOtherwindow);
%                         end
%                     end
                    p=0;
                otherwise
                    error('Unrecognized statistical test type in change.stats.test_type');
            end
            cellByCellStat(k)=p;
        end
        withinCellStats{i,j}=cellByCellStat;
        
        % Output for figure
        switch cset.display_type
            case 'pval x amp'
                withinCellPlotOutput{i,j}=(withinCellStats{i,j}<cset.sigval).*withinCellChanges{i,j};
            case 'amp'
                withinCellPlotOutput{i,j}=withinCellChanges{i,j};
            otherwise
                disp('Using default output: 1 if change during timewindow is significant; 0 otherwise.');
                withinCellPlotOutput{i,j}=withinCellStats{i,j}<cset.sigval;
        end
    end
end

% Get figures for each profile type
rowsOfFig=cell(1,size(withinCellPlotOutput,1));
alltogether=[];
for i=1:size(withinCellPlotOutput,1)
    rowsOfFig{i}=plotCellsByOpto(withinCellPlotOutput,i,1);
    alltogether=[alltogether; [rowsOfFig{i} zeros(size(rowsOfFig{i},1),1)]];
end
alltogether=[alltogether; zeros(1,size(rowsOfFig{i},2)+1)];

% Sort units for display
doSort=0;
if bset.show_profiles(sorting.by_this_behavior)~=1
    disp('Please set the behavioral profile used to sort units to display ON');
    disp('Default: not sorting units');
elseif oset.show_profiles(sorting.by_this_opto)~=1
    disp('Please set the optogenetic profile used to sort units to display ON');
    disp('Default: not sorting units');
else
    doSort=1;
end
if doSort==1
    behNum=sum(bset.show_profiles(1:sorting.by_this_behavior)==1);
    optoNum=sum(oset.show_profiles(1:sorting.by_this_opto)==1);
    sortRow=(behNum-1)*length(optProfiles)+optoNum;
    if strcmp(sorting.order,'descend')
        [~,si]=sort(alltogether(sortRow,1:end-1),2,'descend');
    else
        [~,si]=sort(alltogether(sortRow,1:end-1),2,'ascend');
    end
    alltogether(:,1:end-1)=alltogether(:,si);
    
    % Sort units before returning their activity
    for i=1:size(withinCellResponses,1)
        for j=1:size(withinCellResponses,2)
            temp=withinCellResponses{i,j};
            withinCellResponses{i,j}=temp(si);
            temp=withinCellStats{i,j};
            withinCellStats{i,j}=temp(si);
        end
    end
end
figure();
pcolor(alltogether);
colorbar;
xlabel('Cells');
ylabel('Opto Stim and Behavior');

out.display_matrix=alltogether(1:end-1,1:end-1);
[out.behavior,out.optogenetics,out.change,out.sorting,out.dist,out.traces,out.response]=analysisSettings();

if saveData==1
    d=mkdir([loadDir '\optoTriggeredAnalysis\']);
    if d==1
        save([loadDir '\optoTriggeredAnalysis\withinCellResponses.mat'],'withinCellResponses');
        save([loadDir '\optoTriggeredAnalysis\withinCellStats.mat'],'withinCellStats');
        save([loadDir '\optoTriggeredAnalysis\withinCellAverages.mat'],'withinCellAverages');
        save([loadDir '\optoTriggeredAnalysis\times.mat'],'times');
        save([loadDir '\optoTriggeredAnalysis\optoForProfile.mat'],'optoForProfile');
        save([loadDir '\optoTriggeredAnalysis\out.mat'],'out');
        save([loadDir '\optoTriggeredAnalysis\behForProfile.mat'],'behForProfile');
    end
end

end

function matrix_data=plotCellsByOpto(data,row,suppressFig)

matrix_data=zeros(size(data,2),length(data{row,1}));
for i=1:size(data,2)
    matrix_data(i,:)=data{row,i}';
end

if suppressFig~=1
    figure();
    imagesc(matrix_data);
    colorbar;
    xlabel('Cells');
    ylabel('Opto Stim');
end

end

function outResponses=takeTrialsForEachCell(responses,useTrials)

outResponses=cell(1,length(responses));
for i=1:length(responses)
    r=responses{i};
    outR=r(useTrials,:);
    outResponses{i}=outR;
end
end

function profile=getBehaviorProfile(behavior,opto_stim,b,i)

% b:            b is cell array specifying how to select behavioral profile
% behavior:     array containing some metric of behavior where rows are
%               trials and columns are different time points
% opto_stim:    vector containing the location of the opto stim

% Get behavioral profile specified by b
currb=b{i};
if ~iscell(currb) && ischar(currb)
    % Combines two other behavioral profiles
    % Recursively get these other behavioral profiles, then combine
    startIndex=regexp(currb,'profile');
    firstProfile=str2num(currb(startIndex(1)+7));
    secondProfile=str2num(currb(startIndex(2)+7));
    profile1=getBehaviorProfile(behavior,opto_stim,b{firstProfile});
    profile2=getBehaviorProfile(behavior,opto_stim,b{secondProfile});
    % Combine these profiles
    if ~isempty(regexp(currb,'&','once')) && isempty(regexp(currb,'|','once'))
        % AND
        profile=profile1 & profile2;
    elseif ~isempty(regexp(currb,'|','once')) && isempty(regexp(currb,'&','once'))
        % OR
        profile=profile1 | profile2;
    else
        error(['behavior.profile{' num2str(i) '} should contain & or | but not both']);
    end
elseif ~iscell(currb)
    error(['Incorrect type or format of behavior.profile{' num2str(i) '}']);
elseif iscell(currb)
    if length(currb)==1
        profile=behaviorProfile(behavior,opto_stim,currb{1});
    else
        countSpecs=1;
        countBool=1;
        for j=1:length(currb)
            if iscell(currb{j})
                behSpecs{countSpecs}=currb{j};
                countSpecs=countSpecs+1;
            elseif ischar(currb{j})
                behBool{countBool}=currb{j};
                countBool=countBool+1;
            else
                error(['Incorrect type or format of behavior.profile{' num2str(i) '}']);
            end
        end
        outProfiles=cell(1,length(behSpecs));
        for j=1:length(behSpecs)
            outProfiles{j}=behaviorProfile(behavior,opto_stim,behSpecs{j});
        end
        for j=2:length(behSpecs)
            currBehBool=behBool{j-1};
            if strcmp(currBehBool,'&')
                % AND
                if j==2
                    profile=outProfiles{j-1} & outProfiles{j};
                else
                    profile=profile & outProfiles{j};
                end
            elseif strcmp(currBehBool,'|')
                % OR
                if j==2
                    profile=outProfiles{j-1} | outProfiles{j};
                else
                    profile=profile | outProfiles{j};
                end
            else
                error(['Incorrect type or format of behavior.profile{' num2str(i) '}']);
            end
        end
    end
end
 
end

function profile=behaviorProfile(behavior,opto_stim,b)

behSettings=analysisSettings();

switch behSettings.profile_type
    case 'runningBeforeAndAfterOpto'
        profile=runningProfile(behavior,opto_stim,b);
    otherwise
        error('Unrecognized behavior.profile_type');
end

end

function profile=runningProfile(behavior,opto_stim,b)

% Pass in behavior data as an array where rows are different trials and
% columns are different time points in trial
% Pass in opto_stim structure associated with these trials as a vector
% Low means opto off -- high means opto just finished 
% Note that time points when opto is on have been removed because shutter
% was closed

% Set parameters
runningBefore=b{1}; % Set this to 1 if you want to select trials where animal was running before opto stim 
% Set runningBefore to [0 1] if don't care whether animal was running prior
% to opto stim
runningAfter=b{2}; % Set this to 1 if you want to select trials where animal was stationary prior to opto stim
% Set runningAfter to [0 1] if don't care whether animal was running after
% opto stim
opto_stim_thresh=0.5; % Threshold above which opto stim on 
running_thresh=0.01; % Threshold for behavior above which animal is considered to be running
matchAfterWindowToBeforeWindow=1; % if 1, will match duration of window after opto stim to consider 
% whether mouse was running to duration of window before opto stim;
% otherwise take rest of trial after opto stim

% Get trials associated with this behavioral profile

% Find time points before opto stim
firstOptoStim=find(opto_stim>opto_stim_thresh,1,'first');
beforeInds=1:firstOptoStim-1;

% Find time points after opto stim
if matchAfterWindowToBeforeWindow==1
    if firstOptoStim+length(beforeInds)-1>size(behavior,2)
        afterInds=firstOptoStim:size(behavior,2);
    else
        afterInds=firstOptoStim:firstOptoStim+length(beforeInds)-1;
    end 
else
    afterInds=firstOptoStim:length(opto_stim);
end

% Find trials where animal was running before opto stim
wasRunningBefore=any(abs(behavior(:,beforeInds))>running_thresh,2);
% Find trials where animal was running after opto stim
wasRunningAfter=any(abs(behavior(:,afterInds))>running_thresh,2);
 
% Return trials that fit behavioral profile specified
profile=ismember(wasRunningBefore,runningBefore) & ismember(wasRunningAfter,runningAfter);
end

function profile=optoStimProfile(optoStimTypes,o)

profile=ismember(optoStimTypes,o);

end
