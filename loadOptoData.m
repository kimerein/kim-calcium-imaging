function [opto_shutterTimesRemoved,optoMapping,groups]=loadOptoData(obj,nameOptoCommand,shutterData)
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

% Create a mapping that specifies which opto stim types were given in this
% experiment (i.e., are present in the current data set), and map these
% opto stim types to integer values
[optoMapping,groups]=createOptoMapping(optoData);

% Remove from opto data time points when imaging shutter is closed 
% Note that this will remove voltage commands of opto stimuli from the
% data, because imaging is shuttered while opto is on
% Thus, add delta functions to represent opto stimuli
opto_shutterTimesRemoved=cell(1,size(optoData,1)); % Accomodates different time points shuttered in each trial
for i=1:size(optoData,1)
    [shuttered_opto,shuttered_times]=removeShutteredTimes(obj,shutterData(i,:),optoData(i,:),times,1,30,groups(i));
    % Cut off or fill end of each trial so duration of opto matches duration of
    % imaging acquisition for trial
    opto_shutterTimesRemoved{i}=fixTrialToMatchImaging(obj,shuttered_opto,shuttered_times);
end

end

function physData=fixTrialToMatchImaging(obj,physData,physData_times)

frameDuration=(obj.sabaMetadata.acq.msPerLine/1000)*obj.sabaMetadata.acq.linesPerFrame;
movieDuration=frameDuration*obj.sabaMetadata.acq.numberOfFrames;
physTimeStep=min([mode(diff(physData_times)) physData_times(2)-physData_times(1)]);
physDuration=max(physData_times)+physTimeStep;
if movieDuration<physDuration
    physData=physData(physData_times<=movieDuration);
elseif movieDuration>physDuration
    nElementsToAdd=floor((movieDuration-physDuration)/physTimeStep);
    physData=[physData ones(1,nElementsToAdd).*physData(end)];
end

end

function physData=findPhysData(obj,movieOrder,nameCommand)

% Look for phys command in same path as selected movie files
movName=obj.Movies{1};
dirBreaks=regexp(movName,'\','start');
physPath=movName(1:dirBreaks(end));
listing=dir([physPath nameCommand '*.mat']);
listingnames=cell(1,length(listing));
for i=1:length(listing)
    listingnames{i}=listing(i).name;
end

physData=[];
if ~isempty(listing)
    % Phys command saved during acquisition
    % For each movie file, get associated phys command
    
    for i=1:length(movieOrder)
        ci=movieOrder(i);
        s=obj.Movies{ci};
        parts=regexp(s,'\','split');
        movName=parts{end};
        fi=regexp(movName,'.tif');
        movName=movName(1:fi-1);
        movNumber=movName(end-2:end);
        % Check that expected file is in directory
        physFile=[nameCommand '_' num2str(str2num(movNumber)) '.mat'];
        if any(strcmp(listingnames,physFile))
            warning('off','MATLAB:unknownObjectNowStruct'); % import wave as struct
            % Load phys command
            w=load([physPath physFile]);
            warning('on','MATLAB:unknownObjectNowStruct'); 
            f=fieldnames(w);
            w=w.(f{1});
            if isempty(physData)
                physData=nan(length(movieOrder),length(w.data));
            end
            physData(i,:)=w.data;
        else
            disp(['Missing ' nameCommand ' file from movie directory']);
        end
    end
end

if isempty(physData)
    % No phys data 
    disp('No phys data found');
    times=[];
end
   
end

function [otherData,times]=removeShutteredTimes(obj,shutterData,otherData,times,insertOptoDelta,opto_on_thresh,group)
% Removes shuttered time points from otherData so as to match otherData
% to movies with shuttered frames removed
% 
% times is vector of time points that refers to both shutterData and
% otherData (i.e., shutterData and otherData must have matching lengths 
% and must refer to the same time points)
% shutterData is a vector of values of command voltage to shutter
% otherData is a vector of values
% If insertOptoDelta is true (logical 1), will insert a delta function
% after each shuttered time window IF opto stim was higher than
% opto_on_thresh during that shuttered time window -- these delta functions
% will be added to output otherData

% Get duration of each movie frame
frameDuration=(obj.sabaMetadata.acq.msPerLine/1000)*obj.sabaMetadata.acq.linesPerFrame;

% Find windows when shutter is closed
% Shutter command is assumed to be TTL
shutterStateChanges=find(abs(diff(shutterData))>2.5);
shutterOffWindows=[];
if obj.sabaMetadata.highMeansShutterOff==1
    if shutterData(1)<2.5
        % Shutter starts on
        % Thus shutter off windows are
        for i=1:2:length(shutterStateChanges)
            if i+1>length(shutterStateChanges)
                shutterOffWindows=[shutterOffWindows; shutterStateChanges(i) length(shutterStateChanges)];
            else
                shutterOffWindows=[shutterOffWindows; shutterStateChanges(i) shutterStateChanges(i+1)];
            end
        end
    else
        % Shutter starts off
        % Thus shutter off windows are
        for i=1:2:length(shutterStateChanges)
            if i-1<1
                shutterOffWindows=[shutterOffWindows; 0 shutterStateChanges(1)];
            else
                shutterOffWindows=[shutterOffWindows; shutterStateChanges(i-1) shutterStateChanges(i)];
            end
        end
    end
elseif obj.sabaMetadata.highMeansShutterOff==0
    if shutterData(1)>2.5
        % Shutter starts on
        % Thus shutter off windows are
        for i=1:2:length(shutterStateChanges)
            if i+1>length(shutterStateChanges)
                shutterOffWindows=[shutterOffWindows; shutterStateChanges(i) length(shutterStateChanges)];
            else
                shutterOffWindows=[shutterOffWindows; shutterStateChanges(i) shutterStateChanges(i+1)];
            end
        end
    else
        % Shutter starts off
        % Thus shutter off windows are
        for i=1:2:length(shutterStateChanges)
            if i-1<1
                shutterOffWindows=[shutterOffWindows; 0 shutterStateChanges(1)];
            else
                shutterOffWindows=[shutterOffWindows; shutterStateChanges(i-1) shutterStateChanges(i)];
            end
        end
    end
end
shutterOffWindows=shutterOffWindows+1; % Shift for indexing into times
shutterOffTimes=reshape(times(shutterOffWindows(1:end)),size(shutterOffWindows,1),size(shutterOffWindows,2));

% Take into account mechanical delay in shutter opening
shutterOffTimes(:,2)=shutterOffTimes(:,2)+(obj.sabaMetadata.shutterOpeningTime/1000); % shutterOpeningTime was given in ms, convert to s
% Take into account duration of acquisition of each frame
shutterOffTimes(:,1)=shutterOffTimes(:,1)-frameDuration;

% Remove shuttered time points from otherData
isShutteredTime=zeros(size(times));
for i=1:size(shutterOffTimes,1)
    isShutteredTime(times>=shutterOffTimes(i,1) & times<=shutterOffTimes(i,2))=1;
    if insertOptoDelta==1
        % If otherData is high during this shuttered time window
        if any(otherData(times>=shutterOffTimes(i,1) & times<=shutterOffTimes(i,2))>opto_on_thresh)
            timepointAfterShutter=find(times>shutterOffTimes(i,2),1,'first');
%             % Make first time point after shuttered window high
%             otherData(timepointAfterShutter)=max(otherData(times>=shutterOffTimes(i,1) & times<=shutterOffTimes(i,2)));
            % Make first time point after shuttered window represent which opto
            % stim was presented
            otherData(timepointAfterShutter)=group;
            otherData([1:timepointAfterShutter-1 timepointAfterShutter+1:end])=0;
        end
        
%         % If otherData is low during this shuttered time window (i.e.,
%         % shutter only, no opto stim)
%         if all(otherData(times>=shutterOffTimes(i,1) & times<=shutterOffTimes(i,2))<opto_on_thresh)
%             timepointAfterShutter=find(times>shutterOffTimes(i,2),1,'first');
%             % Make first time point after shuttered window high
%             otherData(timepointAfterShutter)=max(otherData);
%         end
    end
end
otherData=otherData(isShutteredTime==0);
times=times(isShutteredTime==0);
  
end

function [optoMapping,groups]=createOptoMapping(optoData)

% optoData is opto voltage command
% rows are trials; columns are time points

% Get the set of opto stim presented in this expt (i.e., current data set)
noiseThresh=1; % threshold for high-frequency noise
isSorted=zeros(size(optoData,1),1);
groups=nan(size(optoData,1),1);
currGroup=1;
safetyCounter=1;
while any(isSorted==0)
   currInd=find(isSorted==0,1,'first');
   currOpto=optoData(currInd,:);
   groups(currInd)=currGroup;
   for i=1:size(optoData,1)
       if isSorted(i)==0
           if all(abs(currOpto-optoData(i,:))<noiseThresh) % if vectors are equal
               groups(i)=currGroup;
               isSorted(i)=1;
           end
       end
   end
   currGroup=currGroup+1;
   safetyCounter=safetyCounter+1;
   if safetyCounter>size(optoData,1)
       break
   end          
end

% Make a mapping of these opto stim types onto an integer set
uniqueGroups=unique(groups);
optoMapping=cell(length(uniqueGroups),2);
for i=1:length(uniqueGroups)
    % Put in the opto stim type as a vector of values containing voltage
    % command to opto
    f=find(groups==uniqueGroups(i),1,'first');
    optoMapping{i,1}=optoData(f,:);
    % Put in the integer value that maps to this opto stim type
    optoMapping{i,2}=uniqueGroups(i);
end

end

