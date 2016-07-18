function profile=chooseOptoProfile(optoMapping,optoStimTypes,obj)

% Eventually make this a GUI

% For now, just specify manually here
i=1; 
typeNum=optoMapping{i,2};
profile=ismember(optoStimTypes,typeNum);
disp('number of trials in average');
display(sum(profile)); % number of trials in displayed averages

% Plot this opto stim
figure();
realOpto=optoMapping{i,1};
samplingRate=obj.sabaMetadata.phys.settings.inputRate; % Get sampling rate of opto data
times=0:1/samplingRate:(1/samplingRate)*length(realOpto)-(1/samplingRate);
plot(times,realOpto);
title('Opto Stim Type Selected for Average');