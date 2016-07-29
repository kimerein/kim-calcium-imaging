function analyzeTrialByTrial(responses,stats,times)

% Get analysis settings
[bset,oset,cset,~,dist]=analysisSettings();
behProfiles=bset.profiles(bset.show_profiles==1);
optProfiles=oset.profiles(oset.show_profiles==1);

behNum=sum(bset.show_profiles(1:dist.by_this_behavior)==1);
optoNum=sum(oset.show_profiles(1:dist.by_this_opto)==1);

% Get cell responses matching the specified behavioral and optogenetic
% conditions
data=responses{behNum,optoNum};
stat=stats{behNum,optoNum};

% Plot distribution of effects across trials for each cell with a
% significant change after opto stim
data=data(stat<0.05); % Take cells with a significant change
histgrams.x=cell(1,length(data));
histgrams.y=cell(1,length(data));
changesAcrossCells=cell(1,length(data)); 
allchanges=[];
if dist.sortTrials.yes==1
    valAcrossCells=cell(1,length(data)); 
end
for i=1:length(data)
    currData=data{i};
    % Get changes across trials
    if cset.baseline_subtract==1
        changes=nanmean(currData(:,times>=cset.timewindow(1) & times<=cset.timewindow(2)),2)-nanmean(currData(:,times>=cset.baselinewindow(1) & times<=cset.baselinewindow(2)),2);
    else
        changes=nanmean(currData(:,times>=cset.timewindow(1) & times<=cset.timewindow(2)),2);
    end
    if dist.sortTrials.yes==1
        valAcrossCells{i}=nanmean(currData(:,times>=dist.sortTrials.by_window(1) & times<=dist.sortTrials.by_window(2)),2);
    end
    changesAcrossCells{i}=changes;
    allchanges=[allchanges; changes];
end

if dist.sameBins==1
    mi=min(allchanges);
    ma=max(allchanges);
    edges=linspace(mi,ma,dist.nBins+1);
end
for i=1:length(changesAcrossCells)
    changes=changesAcrossCells{i};
    if dist.sameBins==1
        [n,x]=hist(changes,edges);
    else
        [n,x]=hist(changes,dist.nBins);
    end
    if dist.normalize==1
        areaUnder=sum(n);
        n=n./areaUnder; % Normalize so that integral of histogram is 1
    end
    histgrams.x{i}=x;
    histgrams.y{i}=n;
end
% Plot histograms
[h,yOffsets]=plotCurvesOffset([],histgrams.x,histgrams.y,0.1,'k',[]);
title('Distribution of Effects Across Trials -- Each Curve is One Cell');
xlabel('delta F over F');
ylabel('Count');

% Plot effects combining cells
[n,x]=hist(allchanges,dist.nBins);
figure(); 
plot(x,n,'Color','k');
title('Distribution of Effects Across Trials -- Combining All Cells with Sig Effects');
xlabel('delta F over F');
ylabel('Count');

% Further sort trials
if dist.sortTrials.yes==1
    histgrams1.x=cell(1,length(valAcrossCells));
    histgrams2.y=cell(1,length(valAcrossCells));
    for i=1:length(valAcrossCells)
        vals=valAcrossCells{i};
        changes=changesAcrossCells{i};
        switch dist.sortTrials.divide_at
            case 'median'
                changes1=changes(vals>nanmedian(vals));
                changes2=changes(vals<=nanmedian(vals));
            otherwise
                % Use median as default
                changes1=changes(vals>nanmedian(vals));
                changes2=changes(vals<=nanmedian(vals));
        end
        if dist.sameBins==1
            [n1,x1]=hist(changes1,edges);
            [n2,x2]=hist(changes2,edges);
        else
            [n1,x1]=hist(changes1,dist.nBins);
            [n2,x2]=hist(changes,dist.nBins);
        end
        if dist.normalize==1
            areaUnder=sum(n1);
            n1=n1./areaUnder; % Normalize so that integral of histogram is 1
            areaUnder=sum(n2);
            n2=n1./areaUnder; % Normalize so that integral of histogram is 1
        end
        histgrams1.x{i}=x1;
        histgrams1.y{i}=n1;
        histgrams2.x{i}=x2;
        histgrams2.y{i}=n2;
    end
    % Plot overlaid histograms
    h=figure();
    [~,yOffsets]=plotCurvesOffset(h,histgrams1.x,histgrams1.y,0.1,'b',[]);
    plotCurvesOffset(h,histgrams2.x,histgrams2.y,0.1,'r',yOffsets);
    title('Distribution of Effects Across Trials -- Trials Divided Up');
    xlabel('delta F over F');
    ylabel('Count');
end
