function [velocity,continuous_position]=angular_velocity(positionVector)

velThresh=0.03; % max rate reasonable for mouse to run

% Find where wheel angle flips back to 0
% For all times when wheel angle flips back to 0, make curve continuous
continuous_position=positionVector;
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
    