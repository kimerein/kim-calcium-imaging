function plotCaResponse(withinCellAverages,stats,times,optoForProfile)

% Get analysis settings
[bset,oset,cset,~,dist,traces]=analysisSettings();
behProfiles=bset.profiles(bset.show_profiles==1);
optProfiles=oset.profiles(oset.show_profiles==1);

behNum=sum(bset.show_profiles(1:dist.by_this_behavior)==1);
optoNum=sum(oset.show_profiles(1:dist.by_this_opto)==1);

% Get cell responses matching the specified behavioral and optogenetic
% conditions
data=withinCellAverages{behNum,optoNum};
stat=stats{behNum,optoNum};
opto=optoForProfile{behNum,optoNum};

responseAcrossAll=nan(size(data,1),size(data,2));
responseForIncrease=nan(size(data,1),size(data,2));
responseForDecrease=nan(size(data,1),size(data,2));
responseIncludingNonsig=nan(size(data,1),size(data,2));
for i=1:size(data,1)
    currData=data(i,:);
    % Get change
    if cset.baseline_subtract==1
        change=nanmean(currData(times>=cset.timewindow(1) & times<=cset.timewindow(2)),2)-nanmean(currData(times>=cset.baselinewindow(1) & times<=cset.baselinewindow(2)),2);
    else
        change=nanmean(currData(times>=cset.timewindow(1) & times<=cset.timewindow(2)),2);
    end
    if traces.align_to_baseline==1
        currData=currData-nanmean(currData(times>=cset.baselinewindow(1) & times<=cset.baselinewindow(2)),2);
    end
    responseIncludingNonsig(i,:)=currData;
    if stat(i)<0.05
        responseAcrossAll(i,:)=currData;
        if change>0
            responseForIncrease(i,:)=currData;
        elseif change<0
            responseForDecrease(i,:)=currData;
        end
    end
end

optoToPlot=nanmean(opto,1);
optoToPlot=optoToPlot./max(optoToPlot);
% Plot Ca2+ traces across cells
figure();
hax=axes();
plotWStderr(hax,times,responseIncludingNonsig,'k');
xlabel('Time (s)');
ylabel('delta F over F');
title('Response across all cells');
xlim([traces.xlimits(1) traces.xlimits(2)]);
plot(hax,times,optoToPlot.*max(nanmean(responseIncludingNonsig,1)),'Color','c');

figure();
hax=axes();
plotWStderr(hax,times,responseAcrossAll,'k');
xlabel('Time (s)');
ylabel('delta F over F');
title('Response across all cells with a significant change');
xlim([traces.xlimits(1) traces.xlimits(2)]);
plot(hax,times,optoToPlot.*max(nanmean(responseAcrossAll,1)),'Color','c');

figure();
hax=axes();
plotWStderr(hax,times,responseForIncrease,'r');
xlabel('Time (s)');
ylabel('delta F over F');
title('Response across all cells with a significant increase');
xlim([traces.xlimits(1) traces.xlimits(2)]);
plot(hax,times,optoToPlot.*max(nanmean(responseForIncrease,1)),'Color','c');

figure();
hax=axes();
plotWStderr(hax,times,responseForDecrease,'b');
xlabel('Time (s)');
ylabel('delta F over F');
title('Response across all cells with a significant decrease');
xlim([traces.xlimits(1) traces.xlimits(2)]);
plot(hax,times,optoToPlot.*max(nanmean(responseForDecrease,1)),'Color','c');