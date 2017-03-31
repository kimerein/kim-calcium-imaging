function si=plotZscoredData(data,opto,ind,times,tit,beh,si)

% temp=data{ind};
temp=data;

% Transform to Z score
temp=(temp-repmat(nanmean(temp,2),1,size(temp,2)))./repmat(nanstd(temp,[],2),1,size(temp,2));

% Transform to Z score wrt baseline
% baseline=[2.1 3.9];
% temp=(temp-repmat(nanmean(temp(:,times>=baseline(1) & times<=baseline(2)),2),1,size(temp,2)))./repmat(nanstd(temp(:,times>=baseline(1) & times<=baseline(2)),[],2),1,size(temp,2));

% Exclude outliers
% outlierThresh=10^2;
% outlierInd=max(temp,[],2)>outlierThresh;
% temp=temp(outlierInd~=1,:);

% Fill in nans
temp(isnan(temp))=nanmean(nanmean(temp,1),2);

if isempty(si)
    idx=kmeans(temp(:,times>=3 & times<=6),3);
    [~,si]=sort(idx);
    temp=temp(si,:);
    beh=beh(si,:);
else
    temp=temp(si,:);
    beh=beh(si,:);
end

figure(); 
imagesc(temp);
title(tit);

figure();
imagesc(beh);
title('Behavior');

% plot(times,nanmean(temp,1),'Color','b');
% hold on; 
% plot(times,nanmean(temp,1)+nanstd(temp,[],1)./sqrt(size(temp,1)),'Color','b');
% disp('ncells');
% disp(size(temp,1));
% plot(times,nanmean(temp,1)-nanstd(temp,[],1)./sqrt(size(temp,1)),'Color','b');
% % change=nanmean(temp(:,times>=measurewindow(1) & times<=measurewindow(2)),2)-nanmean(temp(:,times>=baseline(1) & times<=baseline(2)),2);
% tempav=nanmean(temp,1);
% tempav(isinf(tempav))=0;
% ma=max(tempav);
% plot(times,nanmean(opto{ind},1).*ma,'Color','r');
% % plot(times,nanmean(beh{ind},1).*ma,'Color',[0.8 0.5 0]);
% plot(times,nanmean(beh,1).*ma,'Color',[0.8 0.5 0]);
% tochange=regexp(tit,'_');
% tit(tochange)=' ';
% title(tit);