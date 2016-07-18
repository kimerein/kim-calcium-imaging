function plotWStderr(hax,x,data,c)
% hax are axes
% x is X axis data
% data is Y axis data
% data rows are different trials, columns are different time points
% c is color

plot(hax,x,nanmean(data,1),'Color',c);
hold on;
plot(hax,x,nanmean(data,1)+nanstd(data,[],1)./sqrt(size(data,1)),'Color',c);
plot(hax,x,nanmean(data,1)-nanstd(data,[],1)./sqrt(size(data,1)),'Color',c);
