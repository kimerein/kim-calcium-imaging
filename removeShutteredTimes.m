function [otherData,times,tstamps]=removeShutteredTimes(obj,shutterData,otherData,times,insertOptoDelta,opto_on_thresh,group,timestamps_nocrop,isOpto)
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
% Make all shutter values above TTL threshold TTL high
shutterData(shutterData>2.5)=5;
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

isShutteredTime=zeros(size(times));
[~,~,~,~,~,~,resp]=analysisSettings();
for i=1:size(shutterOffTimes,1)
    isShutteredTime(times>=shutterOffTimes(i,1) & times<=shutterOffTimes(i,2))=1;
    if isOpto==1
        if insertOptoDelta==1
            if resp.shutter_only~=1
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
            else
                % If otherData is low during this shuttered time window (i.e.,
                % shutter only, no opto stim) FAKE OPTO
                if all(otherData(times>=shutterOffTimes(i,1) & times<=shutterOffTimes(i,2))<opto_on_thresh)
                    timepointAfterShutter=find(times>shutterOffTimes(i,2),1,'first');
                    %             % Make first time point after shuttered window high
                    %             otherData(timepointAfterShutter)=max(otherData);
                    % Make first time point after shuttered window represent which opto
                    % stim was presented
                    otherData(timepointAfterShutter)=group;
                    otherData([1:timepointAfterShutter-1 timepointAfterShutter+1:end])=0;
                end
            end
        end
    end
end

% Remove shuttered time points from otherData
otherData=otherData(isShutteredTime==0);
tstamps=timestamps_nocrop(isShutteredTime==0);
times=times(isShutteredTime==0);
  
end