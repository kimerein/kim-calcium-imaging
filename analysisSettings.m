function varargout=analysisSettings()

% Specify which type of response to use
response.type='raw'; % 'inferred' from CNMF or 'raw'
response.shutter_only=0; % if 1, use shutter only as the stimulus (no opto), else use opto

% Set type of behavior to analyze
behavior.profile_type='runningBeforeAndAfterOpto';
% Options are:
% 'runningBeforeAndAfterOpto'   will divide up behavior based on running
%                               vs. non-running state before/after opto

% Which behavioral profiles to consider for analysis
behavior.profiles{1}={{0; 0}};
behavior.profiles{2}={{1; 1}};
behavior.profiles{3}={{0; 1}};
behavior.profiles{4}={{1; 0}};
behavior.profiles{5}={{0; 0},{1; 1},'|'}; % Specify & or | to combine behavioral profiles
% Will apply logical operators to the profiles in the order that they
% appear
% For example, {{0; 0},{1; 1},{2; 2},{[0 1]; 2},'|','&','|'} will give
% (([0 before, 0 after] OR [1 before, 1 after]) AND [2 before, 2 after]) OR
% [0 or 1 before, 2 after]
behavior.profiles{6}={{[0 1]; [0 1]}};
% To further combine behavioral profiles, make behavior.profiles{i} a
% string specifying which other behavioral profiles to combine
% For example, behavior.profiles{i}='profile5 & profile6' will combine
% behavior.profiles{5} with behavior.profiles{6} using AND
% Only two profiles are allowed to be combined in this manner using | or &

% Specify which behavioral profiles to show in analysis 
% Refers to indices of behavior.profiles
% 1 if show; 0 if don't show
behavior.show_profiles=[1 1 1 1 1 1]; 
if length(behavior.show_profiles)~=length(behavior.profiles)
    error('Length of behavior.show_profiles must match length of behavior.profiles');
end

% Which optogenetic stimulus types to consider for analysis
% optogenetics.profiles{1}=[5 6 14 2 1 7 9 10 11 12 8]; % these are indices into optoMapping; combine all these opto stim types
% optogenetics.profiles{2}=[3 13 4]; 
% optogenetics.profiles{3}=[1:14]; 
% optogenetics.profiles{1}=[2 4 6 7 8 11]; % short stims
% optogenetics.profiles{2}=[1 3 5 9 10 12]; % long stims
% optogenetics.profiles{3}=[1:12]; % all stims
optogenetics.profiles{1}=[1]; % short stims
optogenetics.profiles{2}=[2]; % long stims
optogenetics.profiles{3}=[1:2]; % all stims

% Specify which optogenetic stimulus types to show in analysis
% Refers to indices of optogenetics.profiles
% 1 if show; 0 if don't show
optogenetics.show_profiles=[1 1 1];
if length(optogenetics.show_profiles)~=length(optogenetics.profiles)
    error('Length of optogenetics.show_profiles must match length of optogenetics.profiles');
end

% Change in Ca2+ trace over the course of trial
change.timewindow=[0.7 6.8]; % in seconds with respect to start of trial
change.baselinewindow=[0.255 0.7]; % in seconds with respect to start of trial
% change.timewindow=[6.4 12.4]; % in seconds with respect to start of trial
% change.baselinewindow=[5.955 6.4]; % in seconds with respect to start of trial
% change.timewindow=[7.4 13.5]; % in seconds with respect to start of trial
% change.baselinewindow=[6.955 7.4]; % in seconds with respect to start of trial
% change.baselinewindow=[6.3 7]; % in seconds with respect to start of trial
change.baseline_subtract=1; 
% If baseline_subtract is 1, analysis will subtract mean of Ca2+ trace during baselinewindow from mean during timewindow
% else analysis will return just the mean of Ca2+ trace during timewindow
change.show_windows=1; % if 1, will show time windows on example figure, else suppress this figure
change.stats.test_type='signrank'; % type of statistical test to perform
change.stats.vs_baseline=1; % if 1, will compare the value during timewindow to the value during baselinewindow
change.stats.vs_othertimewindow=[0 1.5]; % if vs_baseline is 1, then will compare the value during timewindow to the value during vs_othertimewindow
change.display_type='pval x amp'; 
% Options are:
% 'pval x amp'  displays the amplitude of the change during timewindow if
%               the p-val is < change.sigval
change.sigval=0.05; 

% How to sort units in display
sorting.by_this_behavior=6; % Refers to indices of behavior.profiles
sorting.by_this_opto=3; % Refers to indices of optogenetics.profiles
% Will sort units according to their responses under these behavior and
% opto conditions
sorting.order='ascend'; % options are 'ascend' or 'descend'

% Show effect distributions for this condition
dist.by_this_behavior=6; % Refers to indices of behavior.profiles
dist.by_this_opto=3; % Refers to indices of optogenetics.profiles
dist.nBins=1000; % Number of bins to use for histograms of effect distributions
dist.normalize=0; % if 1, will normalize histograms
dist.sameBins=1; % if 1, will use same bins for all histograms
dist.sortTrials.yes=1; % if 1, will further divide trials according to Ca2+ trace during window
dist.sortTrials.by_window=[0.255 0.7]; % divide trials according to Ca2+ trace during this window
dist.sortTrials.divide_at='median'; % divide trials at this location; if 'median', will separate trials at the median

% Plot Ca2+ traces
traces.align_to_baseline=1; % if 1, will align traces at baseline (i.e., baseline-subtract)
% baseline window is change.baselinewindow
traces.xlimits=[0.255 12];

% Output
varargout{1}=behavior;
varargout{2}=optogenetics;
varargout{3}=change;
varargout{4}=sorting;
varargout{5}=dist;
varargout{6}=traces;
varargout{7}=response;

end

