function [h,yOffsets]=plotComponentHistograms(data,baselineWindowIndex,timeWindowIndex,nBins,zeroZero,h,yOffsets,c,baselineWindowIndex2,timeWindowIndex2)

% Get histograms for each component
histgrams.x=cell(1,length(data)+1);
histgrams.y=cell(1,length(data)+1);
allchanges=[];
allbases=[];
cellchanges=cell(length(data)+1,1);
for i=1:length(data)
    d=data{i};
    changes=nanmean(d(:,timeWindowIndex),2)-nanmean(d(:,baselineWindowIndex),2);
%     changes=nanmean(d(:,baselineWindowIndex),2);
%     changes=nanmean(d(:,timeWindowIndex),2);
    allchanges=[allchanges; changes];
    allbases=[allbases; nanmean(d(:,baselineWindowIndex),2)];
    % To subsort based on baseline
    med=median(nanmean(d(:,baselineWindowIndex),2));
    changes=changes(nanmean(d(:,baselineWindowIndex),2)<med);
    cellchanges{i}=changes;
    [n,x]=hist(changes,nBins);
    histgrams.x{i}=x;
    histgrams.y{i}=n;
    if zeroZero==1
        temp=n;
        [~,mi]=min(abs(x));
        temp(mi)=nan;
        histgrams.y{i}=temp;
    end
end
med=median(allbases);
allchanges=allchanges(allbases<med);
cellchanges{end}=allchanges;
[n,x]=hist(allchanges,nBins*3);
histgrams.x{end}=x;
histgrams.y{end}=n;


% Plot these histograms together
[h,yOffsets]=plotCurvesOffset(h,histgrams.x,histgrams.y,0.1,c,yOffsets);
    


% % Get histograms for each component
% histgrams.x=cell(1,length(data));
% histgrams.y=cell(1,length(data));
% for i=1:length(data)
%     d=data{i};
%     changes=nanmean(d(:,timeWindowIndex),2)-nanmean(d(:,baselineWindowIndex),2);
%     % To subsort based on baseline
%     med=median(nanmean(d(:,baselineWindowIndex),2));
%     changes=changes(nanmean(d(:,baselineWindowIndex),2)>=med);
%     [n,x]=hist(changes,nBins);
%     histgrams.x{i}=x;
%     histgrams.y{i}=n;
%     if zeroZero==1
%         temp=n;
%         [~,mi]=min(abs(x));
%         temp(mi)=nan;
%         histgrams.y{i}=temp;
%     end
% end
% 
% % Plot these histograms together
% [h,yOffsets]=plotCurvesOffset(h,histgrams.x,histgrams.y,0.1,'r',yOffsets);

baselineWindowIndex=baselineWindowIndex2;
timeWindowIndex=timeWindowIndex2;
% Get histograms for each component
histgrams.x=cell(1,length(data)+1);
histgrams.y=cell(1,length(data)+1);
allchanges=[];
allbases=[];
cellchanges2=cell(length(data)+1,1);
for i=1:length(data)
    d=data{i};
%     changes=nanmean(d(:,timeWindowIndex),2)-nanmean(d(:,fliplr(baselineWindowIndex)),2);
    changes=nanmean(d(:,timeWindowIndex),2)-nanmean(d(:,baselineWindowIndex),2);
%     changes=nanmean(d(:,timeWindowIndex),2);
%     changes=nanmean(d(:,baselineWindowIndex),2);
    allchanges=[allchanges; changes];
    allbases=[allbases; nanmean(d(:,baselineWindowIndex),2)];
    % To subsort based on baseline
    med=median(nanmean(d(:,baselineWindowIndex),2));
    changes=changes(nanmean(d(:,baselineWindowIndex),2)<med);
    cellchanges2{i}=changes;
    [n,x]=hist(changes,nBins);
    histgrams.x{i}=x;
    histgrams.y{i}=n;
    if zeroZero==1
        temp=n;
        [~,mi]=min(abs(x));
        temp(mi)=nan;
        histgrams.y{i}=temp;
    end
end
med=median(allbases);
allchanges=allchanges(allbases<med);
cellchanges2{end}=allchanges;
[n,x]=hist(allchanges,nBins*3);
histgrams.x{end}=x;
histgrams.y{end}=n;

% Plot these histograms together
[h,yOffsets]=plotCurvesOffset(h,histgrams.x,histgrams.y,0.1,'r',yOffsets);

for i=1:length(cellchanges)
    [p,h]=ranksum(cellchanges{i},cellchanges2{i});
    disp(p);
end
disp(nanmedian(cellchanges{end}));
disp(nanmedian(cellchanges2{end}));