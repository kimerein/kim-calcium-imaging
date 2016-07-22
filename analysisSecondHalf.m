function out=analysisSecondHalf(loadDir)

% Load in data from first half of analysis (see
% plot_dFoverF_vs_opto_AsFunctionOfBehavior.m)
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

% Iterate through behavioral profiles
behProfiles={{[0]; [0]}; {[1]; [1]}; {[0]; [1]}; {[1]; [0]}; {[0 1]; [0 1]}};
% behProfiles={{[0]; [0 1]}; {[1]; [0 1]}; {[0]; [1]}; {[1]; [0]}; {[0 1]; [0 1]}};

% Iterate through opto stim types -- ordered here by group num
% For example, run groupOptoStimsGUI from command line to see opto stim
% types and associated group nums (at right)
% From preliminary (behavioral) analysis, it looks like the ranking
% of opto stims in terms of behavioral effects is:
% 250 ms, 1.2 V
% 20 ms, 3 V
% 20 ms, 5 V
% 50 ms, 5 V
% 100 ms, 5 V
% 250 ms, 3 V
% 250 ms, 5 V
% 4 Hz, 5 V
% 12 Hz, 5 V
% 50 Hz, 5 V
% 1 s, 5 V
% optProfiles=[5 6 14 2 7 1 9 10 11 12 8 3 13 4];
% optProfiles={5; 6 ;[14 2] ;[1 7] ;[9 10 11] ;12 ;8 ;3 ;[13 4]};
% Nonrunning to nonrunning
% no fx, increase, no fx, no fx, decrease, decrease, no trials, decrease,
% no trials

% optProfiles={[9 10 11] ;12 ;8 ;3 ;[13 4]};
% optProfiles={[1]};
% optProfiles={[9 10 11 12 8 3 13 4]};
% optProfiles={[1:3 5:12 14]};
optProfiles={[5 6 14 2 1 7 9 10 11 12 8]; [3 13 4]};
% optProfiles={[5 6 14 2 1 7 9 10 11 12 8 3 13 4]};

trialByTrialTimes=trialByTrialTimes-repmat(trialByTrialTimes(:,1),1,size(trialByTrialTimes,2));
times=nanmean(trialByTrialTimes,1);
withinCellAverages=cell(length(behProfiles),length(optProfiles));
% timewindow in seconds with respect to start of trial
timewindow=[2 14];
baselinewindow=[0 1.5];
% timewindow=[2 4];
% timewindow=[2.7 9];
withinCellIntegrals=cell(length(behProfiles),length(optProfiles));
withinCellT=cell(length(behProfiles),length(optProfiles));
for i=1:length(behProfiles)
    b=behProfiles{i};
    % Find trials with desired behavior profile
    behOnlyProfile=behaviorProfile(trialByTrialBeh,avOpto,b);
    
    for j=1:length(optProfiles)
%         o=optProfiles(j);
        o=optProfiles{j};
        % Find trials with desired opto stim type
        optoOnlyProfile=optoStimProfile(optoMapping,optoStimTypes,acq_obj,o);
        % Combine behavioral and opto trial types
        profile=(behOnlyProfile==1) & (optoOnlyProfile==1);
        if sum(profile)==0
            cellByCellAv=zeros(length(useComponents),size(optoTriggeredResponses{1},2));
        else
            % Get average response for each cell
            avResponses=takeTrialsForEachCell(optoTriggeredResponses,profile);
            
            % Get average response across all cells
            cellByCellAv=zeros(length(useComponents),size(avResponses{1},2));
            for k=1:length(useComponents)
                cellByCellAv(k,:)=nanmean(avResponses{useComponents(k)},1);
            end 
        end
        withinCellAverages{i,j}=cellByCellAv;
        
        if i==2 && j==1
            h=figure();
            useBaseInds=zeros(size(times));
            useBaseInds(1:10)=1;
            [~,yOffsets]=plotComponentHistograms(avResponses(useComponents),logical(useBaseInds),times>=timewindow(1) & times<=timewindow(2),100,0,h,[],'k');
        end
        
        if i==2 && j==1
%             plotComponentHistograms(avResponses(useComponents),1:10,times>=timewindow(1) & times<=timewindow(2),100,0,h,yOffsets,'r');
            
            figure();
            hax=axes();
            plotWStderr(hax,times,cellByCellAv,'k');
            hold on;
            plot(hax,times,(avOpto./max(avOpto)).*max(nanmean(cellByCellAv,1)),'Color','c');
            title('Average and Std Err across Cells');
        end
        
        % Get integral of response over timewindow
        withinCellIntegrals{i,j}=nanmean(cellByCellAv(:,times>=timewindow(1) & times<=timewindow(2)),2);
        
        % Get p of ranksum
        cellByCellT=zeros(length(useComponents),1);
        for k=1:length(useComponents)
            currResponse=avResponses{useComponents(k)};
            currIntegrals=nanmean(currResponse(:,times>=timewindow(1) & times<=timewindow(2)),2);
            currBases=nanmean(currResponse(:,1:10),2);
            [p,h,stats]=ranksum(currIntegrals,currBases);
            if isnan(p)
                cellByCellT(k)=nan;
            elseif p<0.05
                cellByCellT(k)=1;
            else
                cellByCellT(k)=nan;
            end
            withinCellT{i,j}=nanmean(cellByCellAv(:,times>=timewindow(1) & times<=timewindow(2)),2).*cellByCellT;
        end
        
%         % Get T statistic of each cell's response
%         cellByCellT=zeros(length(useComponents),1);
%         for k=1:length(useComponents)
%             currResponse=avResponses{useComponents(k)};
%             currIntegrals=sum(currResponse(:,times>=timewindow(1) & times<=timewindow(2)),2);
% %             [~,p,ci,stats]=ttest(currIntegrals);
%             [p,h,stats]=signrank(currIntegrals,0,'method','approximate');
%             if isnan(p)
%                 cellByCellT(k)=0;
% %             elseif p<0.05
%             elseif p<1
% %                 cellByCellT(k)=stats.tstat;
%                 cellByCellT(k)=stats.zval;
%             else
%                 cellByCellT(k)=0;
%             end
% % %             cellByCellT(k)=p;
%         end
%         withinCellT{i,j}=cellByCellT;
    end
end

% Plot integrals as a function of opto stim
nonnon=plotCellsByOpto(withinCellT,1);
title('Non-running to non-running');
runrun=plotCellsByOpto(withinCellT,2);
title('Running to running');
nonrunrun=plotCellsByOpto(withinCellT,3);
title('Non-running to running');
runnonrun=plotCellsByOpto(withinCellT,4);
title('Running to non-running');
allall=plotCellsByOpto(withinCellT,5);
title('All to all');

% nonnon=plotCellsByOpto(withinCellIntegrals,1);
% title('Non-running to non-running');
% runrun=plotCellsByOpto(withinCellIntegrals,2);
% title('Running to running');
% nonrunrun=plotCellsByOpto(withinCellIntegrals,3);
% title('Non-running to running');
% runnonrun=plotCellsByOpto(withinCellIntegrals,4);
% title('Running to non-running');
% allall=plotCellsByOpto(withinCellIntegrals,5);
% title('All to all');

figure();
pcolor([[nonnon zeros(size(nonnon,1),1)]; [runrun zeros(size(nonnon,1),1)]; [nonrunrun zeros(size(nonnon,1),1)]; [runnonrun zeros(size(nonnon,1),1)]; [allall zeros(size(nonnon,1),1)]; zeros(1,size(allall,2)+1)]);
colorbar;
xlabel('Cells');
ylabel('Opto Stim');





end

function matrix_data=plotCellsByOpto(data,row)

matrix_data=zeros(size(data,2),length(data{row,1}));
for i=1:size(data,2)
    matrix_data(i,:)=data{row,i}';
end

figure();
imagesc(matrix_data);
colorbar;
xlabel('Cells');
ylabel('Opto Stim');

end

function outResponses=takeTrialsForEachCell(responses,useTrials)

outResponses=cell(1,length(responses));
for i=1:length(responses)
    r=responses{i};
    outR=r(useTrials,:);
    outResponses{i}=outR;
end
end

function profile=behaviorProfile(behavior,opto_stim,b)

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
opto_stim_thresh=0.75; % Threshold above which opto stim on 
running_thresh=0.1; % Threshold for behavior above which animal is considered to be running
matchAfterWindowToBeforeWindow=1; % if 1, will match duration of window after opto stim to consider 
% whether mouse was running to duration of window before opto stim;
% otherwise take rest of trial after opto stim

% Get trials associated with this behavioral profile

% Find time points before opto stim
firstOptoStim=find(opto_stim>opto_stim_thresh,1,'first');
beforeInds=1:firstOptoStim-1;

% Find time points after opto stim
if matchAfterWindowToBeforeWindow==1
    afterInds=firstOptoStim:firstOptoStim+length(beforeInds)-1;
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

function profile=optoStimProfile(optoMapping,optoStimTypes,obj,o)

% Eventually make this a GUI

% For now, just specify manually here
% i=o; 
% typeNum=optoMapping{i,2};
% profile=ismember(optoStimTypes,typeNum);
profile=ismember(optoStimTypes,o);

% % Plot this opto stim
% figure();
% realOpto=optoMapping{i,1};
% samplingRate=obj.sabaMetadata.phys.settings.inputRate; % Get sampling rate of opto data
% times=0:1/samplingRate:(1/samplingRate)*length(realOpto)-(1/samplingRate);
% plot(times,realOpto);
% title('Opto Stim Type Selected for Average');
end