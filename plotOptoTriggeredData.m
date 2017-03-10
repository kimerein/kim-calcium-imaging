function [change,si]=plotOptoTriggeredData(data,opto,ind,times,tit,beh)

baselinesubtract=1;
% baseline=[1.5 3.5];
baseline=[3 4.224];
measurewindow=[4.224 5.224];

temp=data{ind};

if baselinesubtract==1
    temp=temp-repmat(nanmean(temp(:,times>=baseline(1) & times<=baseline(2)),2),1,size(temp,2));
end

% excludeCells=[18];
excludeCells=[];
takeCells=1:size(temp,1);
temp=temp(takeCells(~ismember(takeCells,excludeCells)),:);

temp=temp(~all(temp==0 | isnan(temp),2),:);

change=nanmean(temp(:,times>=measurewindow(1) & times<=measurewindow(2)),2)-nanmean(temp(:,times>=baseline(1) & times<=baseline(2)),2);
% temp=temp(change>0.01,:);

figure(); 
plot(times,nanmean(temp,1),'Color','b');
hold on; 
plot(times,nanmean(temp,1)+nanstd(temp,[],1)./sqrt(size(temp,1)),'Color','b');
disp('ncells');
disp(size(temp,1));
plot(times,nanmean(temp,1)-nanstd(temp,[],1)./sqrt(size(temp,1)),'Color','b');
% change=nanmean(temp(:,times>=measurewindow(1) & times<=measurewindow(2)),2)-nanmean(temp(:,times>=baseline(1) & times<=baseline(2)),2);
tempav=nanmean(temp,1);
tempav(isinf(tempav))=0;
ma=max(tempav);
plot(times,nanmean(opto{ind},1).*ma,'Color','r');
plot(times,nanmean(beh{ind},1).*ma,'Color',[0.8 0.5 0]);
tochange=regexp(tit,'_');
tit(tochange)=' ';
title(tit);

[~,si]=sort(change);