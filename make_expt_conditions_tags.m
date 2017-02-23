function condition_tags=make_expt_conditions_tags(dataDir)

% Load data from dataDir
load([dataDir '\optoStimTypes.mat']);
load([dataDir '\optoTimestamps.mat']);
load([dataDir '\partwayData_moviematched\acq_obj.mat']);
load([dataDir '\partwayData_moviematched\optoMapping.mat']);
load([dataDir '\partwayData_moviematched\optoTriggeredAnalysis\withinCellStats.mat']);
load([dataDir '\partwayData_moviematched\optoTriggeredAnalysis\withinCellResponses.mat']);
load([dataDir '\partwayData_moviematched\optoTriggeredAnalysis\withinCellAverages.mat']);
hi=load([dataDir '\partwayData_moviematched\optoTriggeredAnalysis\times.mat']);
load([dataDir '\partwayData_moviematched\optoTriggeredAnalysis\out.mat']);
load([dataDir '\partwayData_moviematched\optoTriggeredAnalysis\optoForProfile.mat']);
load([dataDir '\partwayData_moviematched\optoTriggeredAnalysis\behForProfile.mat']);

in.behForProfile=behForProfile;
in.optoForProfile=optoForProfile;
in.out=out;
in.times=hi.times;
in.withinCellAverages=withinCellAverages;
in.withinCellResponses=withinCellResponses;
in.withinCellStats=withinCellStats;
in.optoMapping=optoMapping;
in.optoTimestamps=optoTimestamps;
in.optoStimTypes=optoStimTypes;
in.acq_obj=acq_obj;

% Get description of each experimental condition for comparison across
% experiments
all_condition_tags=cell(size(withinCellResponses));
for i=1:size(withinCellResponses,1)
    in.i=i;
    [condition_tags,condition_tags_pointer]=get_condition_tags_GUI(in);
    condition_tags=condition_tags{1};
    takeThese=condition_tags_pointer(condition_tags_pointer~=0);
%     condition_tags=condition_tags(end-size(withinCellAverages,1)+1:end,end-size(withinCellAverages,2)+1:end);
    condition_tags=condition_tags(takeThese);
    for j=1:length(all_condition_tags(i,:))
        all_condition_tags(i,j)=condition_tags(j);
    end
end
condition_tags=all_condition_tags;
    
save([dataDir '\partwayData_moviematched\optoTriggeredAnalysis\condition_tags.mat'],'condition_tags');