function [opto_stim,optoMovieLength]=alignToMovieFrames(opto_shutterTimesRemoved,C_df,acq_obj)

shutterThresh=0.5; % Assumes that saved shutter signal ranges between 0 and 1

movName=acq_obj.Movies{1};
dirBreaks=regexp(movName,'\','start');
shutterPath=movName(1:dirBreaks(end));
if ~isfield(acq_obj.sabaMetadata,'saveShutterDataFolder')
    error('Need shutter data folder');
else
    listing=dir([shutterPath acq_obj.sabaMetadata.saveShutterDataFolder '\shutterTimesInMovie.mat']);
    if ~isempty(listing)
        a=load([shutterPath acq_obj.sabaMetadata.saveShutterDataFolder '\shutterTimesInMovie.mat']);
        shutterTimesInMovie=a.shutterTimesInMovie;
    else
        error('Need shutter times in movie');
    end
end

% For each trial, force length of opto_shutterTimesRemoved to match length
% of movie
optoMovieLength=cell(1,length(opto_shutterTimesRemoved));
opto_stim=[];
for i=1:length(opto_shutterTimesRemoved)
    optoMovie=zeros(1,length(shutterTimesInMovie{i}));
    % Find end of first shutter epoch in movie
    % Match opto to this frame
    s=shutterTimesInMovie{i};
    shutterOn=find(s>shutterThresh);
    putOptoHere=1;
    for j=1:length(shutterOn)
        % Look for time when shutter opens after first epoch
        if shutterOn(j)+1>length(s)
            % Shutter never reopens
        elseif s(shutterOn(j)+1)<shutterThresh
            putOptoHere=shutterOn(j)+1;
            break
        end
    end
    if ~isnan(putOptoHere)
        optoMovie(putOptoHere)=max(opto_shutterTimesRemoved{i});
    end
    % Remove shuttered times
    optoMovieLength{i}=optoMovie(s<shutterThresh);
    opto_stim=[opto_stim optoMovieLength{i}];
end

if size(C_df,2)~=length(opto_stim)
%     error('Length of opto_stim and movies do not match');
    opto_stim=[opto_stim zeros(1,size(C_df,2)-length(opto_stim))];
end