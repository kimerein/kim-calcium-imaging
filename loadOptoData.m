function [opto_shutterTimesRemoved,optoMapping,groups,optoTimestamps]=loadOptoData(obj,nameOptoCommand,shutterData,movieTimestamps)
% obj is Acquisition2P object from Harvey lab motion correction

nMovies = length(obj.Movies);
movieOrder = 1:nMovies;

if isempty(shutterData)
    shutterData=findPhysData(obj,movieOrder,obj.sabaMetadata.nameShutterCommand);
end
samplingRate=obj.sabaMetadata.phys.settings.inputRate; % Get sampling rate of phys data
times=0:1/samplingRate:(1/samplingRate)*size(shutterData,2)-(1/samplingRate);

% Load in opto data
optoData=findPhysData(obj,movieOrder,nameOptoCommand);
% Get opto data timestamps in relation to movie data
[timestamps,startOfNewTrial]=getDataTimestamps(optoData,times,movieTimestamps);
optoTimestamps.timestamps_nocrop=timestamps;
optoTimestamps.startOfNewTrial_nocrop=startOfNewTrial;

% Check for control (e.g., "bluelightcontrol": opto stim is outside of head) in movie files
areControls=checkForControl(obj);

% Create a mapping that specifies which opto stim types were given in this
% experiment (i.e., are present in the current data set), and map these
% opto stim types to integer values
[optoMapping,groups]=createOptoMapping(optoData,times,areControls);

% Remove from opto data time points when imaging shutter is closed 
% Note that this will remove voltage commands of opto stimuli from the
% data, because imaging is shuttered while opto is on
% Thus, add delta functions to represent opto stimuli
opto_shutterTimesRemoved=cell(1,size(optoData,1)); % Accomodates different time points shuttered in each trial
optoTimestamps.timestamps_aftercrop=cell(1,size(optoData,1));
for i=1:size(optoData,1)
    [shuttered_opto,shuttered_times,tstamps]=removeShutteredTimes(obj,shutterData(i,:),optoData(i,:),times,1,30,groups(i),optoTimestamps.timestamps_nocrop(i,:),1);
    % Cut off or fill end of each trial so duration of opto matches duration of
    % imaging acquisition for trial
    [opto_shutterTimesRemoved{i},optoTimestamps.timestamps_aftercrop{i}]=fixTrialToMatchImaging(obj,shuttered_opto,shuttered_times,tstamps);
end

end

function areControls=checkForControl(obj)

% Parse movie names, looking for key word tag
[~,opto]=analysisSettings();
tag=opto.control_tag; % a string indicating something different about the opto stim in this file, will be grouped separately
areControls=zeros(1,length(obj.Movies));
for i=1:length(obj.Movies)
    currMovName=obj.Movies{i};
    currMovName=fliplr(currMovName);
    startInd=regexp(currMovName,'\','once');
    currMovName=fliplr(currMovName(1:startInd-1));
    if ~isempty(strfind(currMovName,tag))
        % Is a control file
        areControls(i)=1;
    end
end
end

function [optoMapping,groups]=createOptoMapping(optoData,times,areControls)

% optoData is opto voltage command
% rows are trials; columns are time points

% areControls indicates any trials tagged in the movie file name as
% different, group separately

% Get the set of opto stim presented in this expt (i.e., current data set)
[~,oo]=analysisSettings();
noiseThresh=oo.noiseThresh; % threshold for high-frequency noise
timeJitter=oo.acceptable_timeJitter_thresh; % threshold for acceptable time jitter in opto pulse onset/offset
isSorted=zeros(size(optoData,1),1);
groups=nan(size(optoData,1),1);
currGroup=1;
safetyCounter=1;
optoData(areControls==1,:)=optoData(areControls==1,:)+noiseThresh*2;
while any(isSorted==0)
   currInd=find(isSorted==0,1,'first');
   currOpto=optoData(currInd,:);
   groups(currInd)=currGroup;
   for i=1:size(optoData,1)
       if isSorted(i)==0
           if all(abs(currOpto-optoData(i,:))<noiseThresh) % if vectors are equal
               groups(i)=currGroup;
               isSorted(i)=1;
           else
               stillSame=1;
               % Short-cut to bypass the following check in else
               if sum(abs(currOpto-optoData(i,:)))>noiseThresh*(length(currOpto)/2)
                   % Vectors are definitely not the same
                   stillSame=0;
               else
                   % Test whether differences are small time differences in
                   % pulse onset/offset, i.e., acceptable jitter
                   timeStep=times(2)-times(1);
                   differTimes=find(abs(currOpto-optoData(i,:))>noiseThresh);
                   % Find duration of these differing stretches
                   checkedTime=zeros(1,length(currOpto));
                   for j=1:length(differTimes)
                       curr_diffTime=differTimes(j);
                       checkedTime(curr_diffTime)=1;
                       if ismember(curr_diffTime+1,differTimes)
                           % Stretch continues
                           continue
                       else
                           % Stretch ends after this
                           checkedTime(curr_diffTime+1)=1;
                           lengthOfStretch=sum(checkedTime);
                           if lengthOfStretch*timeStep>timeJitter
                               stillSame=0;
                               break
                           else
                               checkedTime=zeros(1,length(currOpto)); % reset
                           end
                       end
                   end
               end
               if stillSame==1 % opto stims should be grouped
                    groups(i)=currGroup;
                    isSorted(i)=1;
               else % opto stims should not be grouped
               end
           end
       end
   end
   currGroup=currGroup+1;
   safetyCounter=safetyCounter+1;
   if safetyCounter>size(optoData,1)
       break
   end          
end
isControlGroups=[];
for i=1:length(groups)
    if min(optoData(i,:))>=noiseThresh*2
        isControlGroups=[isControlGroups groups(i)];
    end
end
isControlGroups=unique(isControlGroups);
optoData(areControls==1,:)=optoData(areControls==1,:)-noiseThresh*2;

% Make a mapping of these opto stim types onto an integer set
uniqueGroups=unique(groups);
optoMapping=cell(length(uniqueGroups),3);
for i=1:length(uniqueGroups)
    % Put in the opto stim type as a vector of values containing voltage
    % command to opto
    f=find(groups==uniqueGroups(i),1,'first');
    optoMapping{i,1}=optoData(f,:);
    % Put in the integer value that maps to this opto stim type
    optoMapping{i,2}=uniqueGroups(i);
    % Code to indicate whether this group is a control
    if ismember(uniqueGroups(i),isControlGroups)
        optoMapping{i,3}='control';
    else
        optoMapping{i,3}='test';
    end
end

end

