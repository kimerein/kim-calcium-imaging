function [velocity,continuous_position]=angular_velocity(positionVector,acq_obj)

samplingRate=acq_obj.sabaMetadata.phys.settings.inputRate; % Get sampling rate of phys data
times=0:1/samplingRate:(1/samplingRate)*length(positionVector)-(1/samplingRate);

velThresh=0.0015; % max rate reasonable for mouse to run
% velThresh=0.001; % max rate reasonable for mouse to run

% Find where wheel angle flips back to 0
% For all times when wheel angle flips back to 0, make curve continuous
% Dealing w a periodic noise source
continuous_position=positionVector;
freqOfNoise=1/(0.128*10^4*(times(2)-times(1))); 
fsResamp=0.128*10^4*freqOfNoise;
pResamp=resample(positionVector,fsResamp,1/(times(2)-times(1)));
newTimes=(0:numel(pResamp)-1)/fsResamp;
avgP=sgolayfilt(pResamp,1,0.128*10^4+1);
continuous_position=avgP;
% Filtering and stitching
filtBin=100;
filtPos=filter(ones(1,filtBin)./filtBin,1,continuous_position);
delayOfFilt=floor((filtBin-1)/2);
filtPos=[filtPos(delayOfFilt:end) zeros(1,length(positionVector)-length(filtPos(delayOfFilt:end)))];
filtPos(1:filtBin)=filtPos(filtBin+1);
filtPos(end-filtBin:end)=filtPos(end-filtBin-1);
continuous_position=filtPos;
flipInds=zeros(1,10);
while ~isempty(flipInds)
    rawDiff=diff(continuous_position);
    flipInds=find(abs(rawDiff)>=velThresh);
    for i=1:length(flipInds)
        if flipInds(i)+1>length(continuous_position)
        else
            continuous_position(flipInds(i)+1:end)=continuous_position(flipInds(i)+1:end)+continuous_position(flipInds(i))-continuous_position(flipInds(i)+1);
        end
    end
end

% Smooth curve to reduce high-frequency noise
% Smooth using a bin size that is small with respect to max velocity
% attainable by mouse, e.g., 10000 for Kim's rig
binSize=1000; 
continuous_position=smooth(continuous_position,binSize);

% Get velocity
velocity=diff(continuous_position,1);
% Fix end of velocity
velocity(end-10+1:end)=velocity(end-10+1).*ones(size(velocity(end-10+1:end)));
    