function [matrixCellResponses,matrixBeh]=selectFinalDatafromBehAndOpto(allTrialsResponses,allTrialsBeh,allOptoTypes,timestep,allOptoData,takeTheseOptoConds,cell_from_expt)

window1=[-2 0]; % wrt opto stim
window2=[0 4]; % wrt opto stim

times=0:timestep:(size(allTrialsResponses{1},2)-1)*timestep;
optoAt=times(nanmean(allOptoData,1)>0);
times=times-optoAt;

matrixBeh=nan(length(allTrialsResponses),length(times));
matrixCellResponses=nan(length(allTrialsResponses),length(times));
for i=1:length(allTrialsResponses)
    curr=allTrialsResponses{i};
    currbeh=allTrialsBeh{i};
    optoProfile=ismember(allOptoTypes{i},takeTheseOptoConds{cell_from_expt{i}});
%     behProfile=constantSpeed_behaviorProfile(currbeh,[-0.0005 0.0005],window1,window2,times);
    behProfile=constantSpeed_behaviorProfile(currbeh,[-1000 1000],window1,window2,times);
    responseProfile=largeResponse_responseProfile(curr,0,window1,window2,times);
    matrixCellResponses(i,:)=nanmean(curr(logical(optoProfile==1 & behProfile'==1),:),1);
%     matrixCellResponses(i,:)=nanmean(curr(logical(optoProfile==1 & behProfile'==1 & responseProfile==1),:),1);
    matrixBeh(i,:)=nanmean(currbeh(logical(optoProfile==1 & behProfile'==1),:),1);
%     matrixBeh(i,:)=nanmean(currbeh(logical(optoProfile==1 & behProfile'==1 & responseProfile==1),:),1);
end  

end

function profile=constantSpeed_behaviorProfile(beh,difference_range,window1,window2,times)

profile=[];

for i=1:size(beh,1)
    currbeh=beh(i,:);
    profile(i)=isConstantSpeed(currbeh,difference_range,window1,window2,times);
end

end

function isconstant=isConstantSpeed(beh,difference_range,window1,window2,times)

beh1=nanmean(beh(times>=window1(1) & times<=window1(2)));
beh2=nanmean(beh(times>=window2(1) & times<=window2(2)));
if abs(beh1-beh2)>difference_range(1) && abs(beh1-beh2)<difference_range(2)
    isconstant=1;
else
    isconstant=0;
end

end

function profile=largeResponse_responseProfile(response,threshold,window1,window2,times)

resp1=nanmean(response(:,times>=window1(1) & times<=window1(2)),2);
resp2=nanmean(response(:,times>=window2(1) & times<=window2(2)),2);
profile=(resp2-resp1)>threshold;

end