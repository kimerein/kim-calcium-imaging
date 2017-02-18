    function [opto_shutterTimesRemoved,optoTimestamps]=loadGenericPhysData(obj,nameOptoCommand,shutterData,movieTimestamps)
% obj is Acquisition2P object from Harvey lab motion correction

nMovies = length(obj.Movies);
movieOrder = 1:nMovies;

if isempty(shutterData)
    shutterData=findPhysData(obj,movieOrder,obj.sabaMetadata.nameShutterCommand);
end
samplingRate=obj.sabaMetadata.phys.settings.inputRate; % Get sampling rate of phys data
times=0:1/samplingRate:(1/samplingRate)*size(shutterData,2)-(1/samplingRate);

% Load in phys data
optoData=findPhysData(obj,movieOrder,nameOptoCommand);
% Get phys data timestamps in relation to movie data
[timestamps,startOfNewTrial]=getDataTimestamps(optoData,times,movieTimestamps);
optoTimestamps.timestamps_nocrop=timestamps;
optoTimestamps.startOfNewTrial_nocrop=startOfNewTrial;

% Remove from phys data time points when imaging shutter is closed 
opto_shutterTimesRemoved=cell(1,size(optoData,1)); % Accomodates different time points shuttered in each trial
optoTimestamps.timestamps_aftercrop=cell(1,size(optoData,1));
for i=1:size(optoData,1)
    [shuttered_opto,shuttered_times,tstamps]=removeShutteredTimes(obj,shutterData(i,:),optoData(i,:),times,1,0.5,[],optoTimestamps.timestamps_nocrop(i,:),0);
    % Cut off or fill end of each trial so duration of opto matches duration of
    % imaging acquisition for trial
    [opto_shutterTimesRemoved{i},optoTimestamps.timestamps_aftercrop{i}]=fixTrialToMatchImaging(obj,shuttered_opto,shuttered_times,tstamps);
end

end