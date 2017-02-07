function [opto_stim,optoMovieLength]=alignToMovieFrames(opto_shutterTimesRemoved,C_df,acq_obj,isOpto)

% if isOpto==1, will zero out all values except when opto on (just as shutter
% closes)

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

nMovies = length(acq_obj.Movies);
movieOrder = 1:nMovies;
movieOrder([1 acq_obj.motionRefMovNum]) = [acq_obj.motionRefMovNum 1];

samplingRate=acq_obj.sabaMetadata.phys.settings.inputRate; % Get sampling rate of phys data
times=0:1/samplingRate:(1/samplingRate)*length(opto_shutterTimesRemoved{1})-(1/samplingRate);

% Get times associated with each movie frame
frameDuration=(acq_obj.sabaMetadata.acq.msPerLine/1000)*acq_obj.sabaMetadata.acq.linesPerFrame;
movieDuration=frameDuration*acq_obj.sabaMetadata.acq.numberOfFrames;
movieTimes=0:frameDuration:frameDuration*acq_obj.sabaMetadata.acq.numberOfFrames-frameDuration;

% For each trial, force length of opto_shutterTimesRemoved to match length
% of movie
optoMovieLength=cell(1,length(opto_shutterTimesRemoved));
opto_stim=[];
[~,~,~,~,~,~,resp]=analysisSettings();
for i=1:length(opto_shutterTimesRemoved)
    optoMovie=zeros(1,length(shutterTimesInMovie{i}));
    s=shutterTimesInMovie{i};
    if isOpto==1
        if resp.shutter_only~=1
            % Find end of shutter epoch in movie that is closest in time to
            % beginning of opto pulse
            beginning_optopulse=times(find(opto_shutterTimesRemoved{i}>0.3,1,'first'));
            s=shutterTimesInMovie{i};
            shutterOn=find(s>shutterThresh);
            timeToOptoPulse=movieTimes(shutterOn)-beginning_optopulse;
            [~,closestToOpto]=min(abs(timeToOptoPulse));
            
            % Match opto to this frame
            putOptoHere=1;
            for j=closestToOpto:length(shutterOn)
                % Look for time when shutter opens after first epoch containing
                % opto stim
                if shutterOn(j)+1>length(s)
                    % Shutter never reopens
                elseif s(shutterOn(j)+1)<shutterThresh
                    putOptoHere=shutterOn(j)+1;
                    break
                end
            end
        else
            % Find end of last shutter epoch in movie
            % Put fake opto here
            s=shutterTimesInMovie{i};
            shutterOff=find(s<shutterThresh);
            putOptoHere=length(s);
            for j=length(shutterOff):-1:1
                % Look for time when shutter last closes
                if shutterOff(j)-1<1
                    % Shutter never closed
                elseif s(shutterOff(j)-1)>shutterThresh
                    putOptoHere=shutterOff(j);
                    break
                end
            end
        end
    end
    if isOpto==1
        if ~isnan(putOptoHere)
            optoMovie(putOptoHere)=max(opto_shutterTimesRemoved{i});
        end
    else
        % Force length of opto_shutterTimesRemoved{i} to match length of
        % movie
        resampleRate=floor(length(opto_shutterTimesRemoved{i})/length(shutterTimesInMovie{i}));
        temp=decimate(opto_shutterTimesRemoved{i},resampleRate);
        optoMovie=temp(1:length(shutterTimesInMovie{i}));
    end
    % Remove shuttered times 
    optoMovieLength{i}=optoMovie(s<shutterThresh);
    if length(optoMovieLength{i})~=acq_obj.correctedMovies.slice.channel.size(i,3)
        disp('Mismatch between optoMovie length and length of actual movie');
        disp('in alignToMovieFrames.m');
        % Patch BUT THE ABOVE SHOULD NOT HAPPEN
        if acq_obj.correctedMovies.slice.channel.size(i,3)>length(optoMovieLength{i})
            optoMovieLength{i}=[optoMovieLength{i} zeros(1,acq_obj.correctedMovies.slice.channel.size(i,3)-length(optoMovieLength{i}))];
        elseif acq_obj.correctedMovies.slice.channel.size(i,3)<length(optoMovieLength{i})
            temp=optoMovieLength{i};
            optoMovieLength{i}=temp(1:acq_obj.correctedMovies.slice.channel.size(i,3));
        end
    end
    opto_stim=[opto_stim optoMovieLength{i}];
end