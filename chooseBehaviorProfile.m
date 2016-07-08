function profile=chooseBehaviorProfile(behavior,opto_stim)

% Pass in behavior data as an array where rows are different trials and
% columns are different time points in trial
% Pass in opto_stim structure associated with these trials as a vector
% Low means opto off -- high means opto just finished 
% Note that time points when opto is on have been removed because shutter
% was closed

% Set parameters
runningBefore=[0]; % Set this to 1 if you want to select trials where animal was running before opto stim 
% Set runningBefore to [0 1] if don't care whether animal was running prior
% to opto stim
runningAfter=[0]; % Set this to 1 if you want to select trials where animal was stationary prior to opto stim
% Set runningAfter to [0 1] if don't care whether animal was running after
% opto stim
opto_stim_thresh=0.5; % Threshold above which opto stim on (opto stim should now range between 0 and 1)
running_thresh=0.1; % Threshold for behavior above which animal is considered to be running
matchAfterWindowToBeforeWindow=1; % if 1, will match duration of window after opto stim to consider 
% whether mouse was running to duration of window before opto stim;
% otherwise take rest of trial after opto stim

% Get trials associated with this behavioral profile

% Find time points before opto stim
firstOptoStim=find(opto_stim>opto_stim_thresh,1,'first');
beforeInds=1:firstOptoStim-1;

% Find time points after opto stim
firstOptoStim=find(opto_stim>opto_stim_thresh,1,'last');
if matchAfterWindowToBeforeWindow==1
    afterInds=firstOptoStim:firstOptoStim+length(beforeInds)-1;
else
    afterInds=firstOptoStim:length(opto_stim);
end

% Find trials where animal was running before opto stim
wasRunningBefore=any(abs(behavior(:,beforeInds))>running_thresh,2);
% Find trials where animal was running after opto stim
wasRunningAfter=any(abs(behavior(:,afterInds))>running_thresh,2);

% Return trials that fit behavioral profile specified
profile=ismember(wasRunningBefore,runningBefore) & ismember(wasRunningAfter,runningAfter);