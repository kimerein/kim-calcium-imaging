function [physData,new_tstamps]=fixTrialToMatchImaging(obj,physData,physData_times,tstamps)

frameDuration=(obj.sabaMetadata.acq.msPerLine/1000)*obj.sabaMetadata.acq.linesPerFrame;
movieDuration=frameDuration*obj.sabaMetadata.acq.numberOfFrames;
physTimeStep=min([mode(diff(physData_times)) physData_times(2)-physData_times(1)]);
physDuration=max(physData_times)+physTimeStep;
if movieDuration<physDuration
    physData=physData(physData_times<=movieDuration);
    new_tstamps=tstamps(physData_times<=movieDuration);
elseif movieDuration>physDuration
    nElementsToAdd=floor((movieDuration-physDuration)/physTimeStep);
    physData=[physData ones(1,nElementsToAdd).*physData(end)];
    new_tstamps=[tstamps tstamps(end)+physTimeStep:physTimeStep:tstamps(end)+physTimeStep*nElementsToAdd];
end

end