function [matrixCellResponses,matrixBeh,selectedAverageCellResponses,behSelected,allbehprofiles,allTrialsResponses,allTrialsBeh,alloptoprofiles]=selectDataByBehavior(allSingleTrialBeh,allCellTrialResponses,allOptoData,timestep,optoStimTypes,useThese,cell_from_expt)

window1=[-2 0]; % wrt opto stim
window2=[0 4]; % wrt opto stim

fix=1;

if fix==1
    for i=1:length(allSingleTrialBeh)
        allSingleTrialBeh{i}=fixRunningData(allSingleTrialBeh{i});
    end
end

times=0:timestep:(size(allSingleTrialBeh{1},2)-1)*timestep;
optoAt=times(nanmean(allOptoData,1)>0);
times=times-optoAt;

selectedAverageCellResponses=cell(1,length(allCellTrialResponses));
behSelected=cell(1,length(allCellTrialResponses));
matrixCellResponses=nan(length(allCellTrialResponses),length(times));
matrixBeh=nan(length(allCellTrialResponses),length(times));
allbehprofiles=cell(1,length(allCellTrialResponses));
allTrialsResponses=cell(1,length(allCellTrialResponses));
allTrialsBeh=cell(1,length(allCellTrialResponses));
alloptoprofiles=cell(1,length(allCellTrialResponses));;
for i=1:length(allCellTrialResponses)
    curr=allCellTrialResponses{i};
    currbeh=allSingleTrialBeh{i};
    currsopto=optoStimTypes{cell_from_expt{i}}; 
    optoProfile=ismember(currsopto,useThese{cell_from_expt{i}});
    % Trim behavior profile according to useThese
    alloptoprofiles{i}=currsopto;
    currbeh=currbeh(optoProfile==1,:);
    if size(curr,1)~=size(currbeh,1) && length(optoProfile)==size(curr,1)
        curr=curr(optoProfile==1,:);
    end
    if length(optoProfile)~=size(curr,1)
        alloptoprofiles{i}=currsopto(optoProfile==1);
    end
    if size(curr,1)~=size(currbeh,1)
        disp(['Error in sizing of i = ' num2str(i)]); 
        continue
    end
    profile=constantSpeed_behaviorProfile(currbeh,[-0.0005 0.0005],window1,window2,times);
    selectedAverageCellResponses{i}=nanmean(curr(profile==1,:),1);
    behSelected{i}=nanmean(currbeh(profile==1,:),1);   
    allbehprofiles{i}=profile;
    allTrialsResponses{i}=curr;
    allTrialsBeh{i}=currbeh;
    matrixCellResponses(i,:)=selectedAverageCellResponses{i};
    matrixBeh(i,:)=behSelected{i};
end

figure(); 
plot(times,nanmean(matrixCellResponses,1),'Color','k');
hold on;
plot(times,nanmean(matrixBeh,1),'Color',[0.2 0.2 0.2]);
disp(size(matrixCellResponses,1));

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