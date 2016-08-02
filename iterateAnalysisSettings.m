function [behavior,optogenetics,change,sorting,dist,traces,response]=iterateAnalysisSettings(behavior,optogenetics,change,sorting,dist,traces,response,clearPersistentCounter)

persistent counter

if clearPersistentCounter==1
    counter=1;
    return
end

if isempty(counter)
    counter=1;
end

% Iterate analysis on the following conditions, saving results for each
% Real opto stim
% Fake opto stim (i.e., shutter only)
% 
% For each ... iterate analysis on all behavioral conditions
% By all opto stim conditions

switch counter
    % I. Real opto stim
    case 1
        response.shutter_only=0; % if 1, use shutter only as the stimulus (no opto), else use opto
        change.timewindow=[0.7 6.8]; % in seconds with respect to start of trial
        change.baselinewindow=[0.255 0.7]; % in seconds with respect to start of trial
    % II. Shutter only
    case 2
        response.shutter_only=1; % if 1, use shutter only as the stimulus (no opto), else use opto
        change.timewindow=[7.4 13.5]; % in seconds with respect to start of trial
        change.baselinewindow=[6.955 7.4]; % in seconds with respect to start of trial
    otherwise
        % Have default be real opto stim
        response.shutter_only=0; % if 1, use shutter only as the stimulus (no opto), else use opto
        change.timewindow=[0.7 6.8]; % in seconds with respect to start of trial
        change.baselinewindow=[0.255 0.7]; % in seconds with respect to start of trial
end

counter=counter+1;
    
    
    
