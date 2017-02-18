function [timestamps,startOfNewTrial]=getDataTimestamps(data,trialtimes,movieTimestamps)

% movieTimestamps contains
% movieTimestamps.timestamps
% movieTimestamps.startOfNewTrial
% movieTimestamps.timestamps_nocrop
% movieTimestamps.startOfNewTrial_nocrop

% First test whether number of trials is equivalent between movie and data
% If not, throw error
if size(data,1)~=sum(movieTimestamps.startOfNewTrial_nocrop==1)
    disp('error in getDataTimestamps');
    error('data and movie have different numbers of trials');
end

% From trialtimes, which are the data times with respect to start of each
% trial, get the times of each data point with respect to the movie times
% Align timestamps at the start of each trial
startOfNewTrial=zeros(size(data));
trialStarts=movieTimestamps.timestamps_nocrop(movieTimestamps.startOfNewTrial_nocrop==1);
timestamps=repmat(trialtimes,size(data,1),1)+repmat(reshape(trialStarts,length(trialStarts),1),1,size(data,2));
startOfNewTrial(:,1)=1;