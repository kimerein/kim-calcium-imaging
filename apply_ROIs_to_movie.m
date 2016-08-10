function traces=apply_ROIs_to_movie(ROIs,acq_obj)

% ROIs are a cell array of ROI masks
% The value at each point in a mask indicates the weighting of activity at
% that pixel for determing spatial extent of ROI

% acq_obj contains Acquisition2P object (output of motion correction)
% with pointers to motion-corrected movies

% Load Acquisition2P object (output of motion correction)
a=load(acq_obj);
names=fieldnames(a);
obj=a.(names{1});

% Apply ROIs to motion-corrected movie to get activity traces over time
% Note: no background or neuropil subtraction
castType='single';
sliceNum=1;
channelNum=obj.motionRefChannel;
traces=cell(length(ROIs),length(obj.correctedMovies.slice.channel.fileName)); % number of ROIs by number of trials
for i=1:length(obj.correctedMovies.slice.channel.fileName)
    disp(i);
    mov=readCor(obj,i,castType,sliceNum,channelNum);
    % mov is the motion-corrected movie
    % as a 3D array where the third dimension is time
    % and first 2 dimensions are spatial
    for j=1:length(ROIs)
        currTrace=nan(1,size(mov,3));
        currROI=ROIs{j};
        for k=1:size(mov,3)
            currTrace(k)=sum(sum(currROI.*mov(:,:,k)));
        end
        traces{j,i}=currTrace;
    end
end
        
figure();
plot(traces{1,1});