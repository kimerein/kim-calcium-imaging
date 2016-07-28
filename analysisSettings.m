function varargout=analysisSettings()

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
optogenetics.profiles{1}=[5 6 14 2 1 7 9 10 11 12 8]; % these are indices into optoMapping; combine all these opto stim types
optogenetics.profiles{2}=[3 13 4]; % these are indices into optoMapping; combine all these opto stim types

% Specify which optogenetic stimulus types to show in analysis
% Refers to indices of optogenetics.profiles
% 1 if show; 0 if don't show
optogenetics.show_profiles=[1 1];
if length(optogenetics.show_profiles)~=length(optogenetics.profiles)
    error('Length of optogenetics.show_profiles must match length of optogenetics.profiles');
end

% Change in Ca2+ trace over the course of trial
change.timewindow=[2 6]; % in seconds with respect to start of trial
change.baselinewindow=[0 1.5]; % in seconds with respect to start of trial
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

% Output
varargout{1}=behavior;
varargout{2}=optogenetics;
varargout{3}=change;

end

