function varargout=analysisSettings()

% Note that this analysis currently assumes the following:
% 1. The "fake opto stim" (i.e., shutter without opto stim) is the last
% shuttered epoch in each trial
% 2. The opto stim occurs only once in each trial (although the opto stim
% may itself be a pulse pattern)

% Specify which type of response to use
response.F_isSlidingWindow=true; % if true, F in delta F over F is from a sliding window
response.F_nTrials=20; % number of trials to use in sliding window to calculate F
response.ITI=2.5; % inter-trial interval in seconds
response.discardTheseFiles=[]; % used for debugging, arbitrary daq output value may cause trouble, will discard these file numbers from further analysis
% in terms of index into movies list in Acquisition2P object
response.type='inferred'; % 'inferred' from CNMF, 'raw' Ca2+ traces with only baseline subtraction, or 'spikes' from deconvolution
response.shutter_only=0; % if 1, use shutter only as the stimulus (no opto), else use opto
response.include_previous_trial=0; % if 1, will add some data from end of previous trial to beginning of each trial
response.n_previous_trial_inds=10; % how many indices from previous trial end to add to beginning of each trial

% Any control files should be labeled with a tag in the movie file name
% For example, the tag "bluelightcontrol" might indicate that the
% optogenetic stimulus was positioned outside of the head
optogenetics.control_tag='bluelightcontrol'; % If this string is in the movie file name, this optogenetic stimulus will be treated as
                                             % different from an identical-looking optogenetic stimulus lacking this tag in the movie file name

% Set type of behavior to analyze
behavior.profile_type='runningBeforeAndAfterOpto';
% Options are:
% 'runningBeforeAndAfterOpto'   will divide up behavior based on running
%                               vs. non-running state before/after opto

% Noise thresholds for automatically grouping opto stims - User will then have the option
% to manually change this grouping
optogenetics.noiseThresh=50; % High-frequency noise threshold (below this, difference is just high-freq noise)
optogenetics.acceptable_timeJitter_thresh=0.0009; % in seconds, if opto pulses differ in onset/offset timing by less than this amount, 
                                                  % they will be considered the same stim
                                                 
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
% these numbers are indices into optoMapping; combine all these opto stim types
optoStimTypes=[1:6]; % all opto stim types in this expt
optogenetics.profiles{1}=[1:6]; 
% optogenetics.profiles{2}=[1:8 10]; 
% optogenetics.profiles{3}=[9]; 
% optogenetics.profiles{4}=[1:3]; % all opto stim types
% Also consider each opto stim type individually
startLength=length(optogenetics.profiles)+1;
j=1;
for i=startLength:startLength+length(optoStimTypes)-1
    optogenetics.profiles{i}=optoStimTypes(j);
    j=j+1;
end

% Specify which optogenetic stimulus types to show in analysis
% Refers to indices of optogenetics.profiles
% 1 if show; 0 if don't show
% optogenetics.show_profiles=[1 1 1 1 1 1 1 1];
optogenetics.show_profiles=ones(size(optogenetics.profiles)); 
% optogenetics.show_profiles=ones(1,length(optogenetics.profiles));
if length(optogenetics.show_profiles)~=length(optogenetics.profiles)
    error('Length of optogenetics.show_profiles must match length of optogenetics.profiles');
end

% Change in Ca2+ trace over the course of trial
change.timewindow=[3.9 3.9+6]; % in seconds with respect to start of trial
change.baselinewindow=[3.9-2 3.9]; % in seconds with respect to start of trial
% change.timewindow=[6.4 12.4]; % in seconds with respect to start of trial
% change.baselinewindow=[5.955 6.4]; % in seconds with respect to start of trial
% change.timewindow=[7.4 13.5]; % in seconds with respect to start of trial
% change.baselinewindow=[6.955 7.4]; % in seconds with respect to start of trial
% change.baselinewindow=[6.3 7]; % in seconds with respect to start of trial
change.baseline_subtract=0; 
% If baseline_subtract is 1, analysis will subtract mean of Ca2+ trace during baselinewindow from mean during timewindow
% else analysis will return just the mean of Ca2+ trace during timewindow
change.show_windows=1; % if 1, will show time windows on example figure, else suppress this figure
change.stats.test_type='signrank'; % type of statistical test to perform
change.stats.vs_baseline=1; % if 1, will compare the value during timewindow to the value during baselinewindow
change.stats.vs_othertimewindow=[0 1.5]; % if vs_baseline is 0, then will compare the value during timewindow to the value during vs_othertimewindow
change.display_type='amp'; 
% Options are:
% 'pval x amp'  displays the amplitude of the change during timewindow if
%               the p-val is < change.sigval
% 'amp'         displays the amplitude of the change during timewindow
change.sigval=0.05; 
% change.sigval=1; 

% How to sort units in display
sorting.by_this_behavior=5; % Refers to indices of behavior.profiles
sorting.by_this_opto=1; % Refers to indices of optogenetics.profiles
% Will sort units according to their responses under these behavior and
% opto conditions
sorting.order='ascend'; % options are 'ascend' or 'descend'

% Show effect distributions for this condition
% In plotCaResponse
dist.by_this_behavior=6; % Refers to indices of behavior.profiles
dist.by_this_opto=1; % Refers to indices of optogenetics.profiles
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

% If am currently iterating multiple analyses, allow
% iterateAnalysisSettings.m to modify these settings before returning
iterate.yes=0; % If 1, will allow iterateAnalysisSettings.m to modify settings values

if iterate.yes==1
    [behavior,optogenetics,change,sorting,dist,traces,response]=iterateAnalysisSettings(behavior,optogenetics,change,sorting,dist,traces,response,0,0);
end

% Output
varargout{1}=behavior;
varargout{2}=optogenetics;
varargout{3}=change;
varargout{4}=sorting;
varargout{5}=dist;
varargout{6}=traces;
varargout{7}=response;

end