function selectDataByBehavior(allSingleTrialBeh,allCellTrialResponses,timestep)

times=0:timestep:(size(allSingleTrialBeh{1},2)-1)*timestep;

selectedAverageCellResponses=cell(1,length(allCellTrialResponses));




end

function profile=constantSpeed_behaviorProfile(beh,difference_range,window1,window2,times)

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