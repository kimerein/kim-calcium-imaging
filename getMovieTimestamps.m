function [timestamps,startOfNewTrial,timestamps_nocrop,startOfNewTrial_nocrop,shutterTimesInMovie]=getMovieTimestamps(shutterData_dir,acq_obj,ITI)

% Load in shutterTimesInMovie
a=load([shutterData_dir '\shutterTimesInMovie.mat']);
shutterTimesInMovie=a.shutterTimesInMovie;

howmany_times=0;
howmany_times_uncropped=0;
for i=1:length(shutterTimesInMovie)
    howmany_times=howmany_times+sum(shutterTimesInMovie{i}==0);
    howmany_times_uncropped=howmany_times_uncropped+length(shutterTimesInMovie{i});
end

timestamps=nan(1,howmany_times);
startOfNewTrial=zeros(1,howmany_times);
timestamps_nocrop=nan(1,howmany_times_uncropped);
startOfNewTrial_nocrop=zeros(1,howmany_times_uncropped);
% Find movie times if had not removed some frames from movie
frameDuration=(acq_obj.sabaMetadata.acq.msPerLine/1000)*acq_obj.sabaMetadata.acq.linesPerFrame;
movieDuration=frameDuration*acq_obj.sabaMetadata.acq.numberOfFrames;
movieTimes=0:frameDuration:frameDuration*acq_obj.sabaMetadata.acq.numberOfFrames-frameDuration;
movieStep=movieTimes(2)-movieTimes(1);
j=1;
currentTime=0;
j_nocrop=1;
for i=1:length(shutterTimesInMovie)
    startOfNewTrial(j)=1;
    startOfNewTrial_nocrop(j_nocrop)=1;
    % Now fill in timestamps for this movie
    movieTimesNow_noCrop=movieTimes+currentTime;
    movieTimesNow_afterCrop=movieTimesNow_noCrop(shutterTimesInMovie{i}==0);
    timestamps(j:j+sum(shutterTimesInMovie{i}==0)-1)=movieTimesNow_afterCrop;
    timestamps_nocrop(j_nocrop:j_nocrop+length(shutterTimesInMovie{i})-1)=movieTimesNow_noCrop;
    % Increment index and currentTime
    j=j+sum(shutterTimesInMovie{i}==0);
    j_nocrop=j_nocrop+length(shutterTimesInMovie{i});
    currentTime=movieTimesNow_noCrop(end)+ITI;
end