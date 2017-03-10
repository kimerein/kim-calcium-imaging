function [allCellResponses,allOptoData,condition_labels,returnTimes,allBehData,timestep]=align_cellResponses(alldata,doAverage)

% This assumes that sampling rate is constant across expts -- will fix this
% later

for i=1:length(alldata.condition_tags)
    temp=alldata.condition_tags{i};
    alldata.condition_tags{i}=temp{1};
end
condition_labels=alldata.condition_tags;
allCellResponses=cell(1,length(alldata.optoForProfile));
allOptoData=cell(1,length(alldata.optoForProfile));
allBehData=cell(1,length(alldata.optoForProfile));
if doAverage==1
    % cycle through expt conditions
    alltimes=nan(1,length(alldata.times));
    for i=1:length(alldata.times)
        currtimes=alldata.times{i};
        alltimes(i)=currtimes(2)-currtimes(1);
    end
    for i=1:length(alldata.optoForProfile)
        temp=alldata.optoForProfile{i};
        tempbeh=alldata.behForProfile{i};
        celltemp=alldata.withinCellAverages{i};
        optoExpts=cell(1,length(temp));
        behExpts=cell(1,length(temp));
        exptTimes=cell(1,length(temp));
        cellAvs=cell(1,length(temp));
        nCells=0;
        for j=1:length(temp)
            optoStims=temp{j};
            behStims=tempbeh{j};
            av=nanmean(optoStims,1);
            avbeh=nanmean(behStims,1);
            optoExpts{j}=av;
            behExpts{j}=avbeh;
            exptTimes{j}=alldata.times{j};
            withinExpt_cellAvs=celltemp{j};
            cellAvs{j}=withinExpt_cellAvs;
            nCells=nCells+size(withinExpt_cellAvs,1);
        end
        if length(temp)>1       % data from more than one expt for this condition
            % need to align to opto stim
            optoAt=nan(1,length(temp));
            opto_lengths=nan(1,length(temp));
            for j=1:length(temp)
                optoStims=optoExpts{j};
                [~,ma]=max(optoStims);
                optoAt(j)=ma;
                opto_lengths(j)=length(optoStims);
            end
            % find latest opto stim -- will pad all others to this length
            ma=max(optoAt);
            mi=min(optoAt);
            allOptos=nan(length(temp),ma+(max(opto_lengths)-mi));
            allResponses=nan(nCells,size(allOptos,2));
            allBeh=nan(length(temp),ma+(max(opto_lengths)-mi));
            % align opto stims and cell responses
            onCell=1;
            for j=1:length(temp)
                allOptos(j,ma-optoAt(j)+1:ma-optoAt(j)+opto_lengths(j))=optoExpts{j};
                allBeh(j,ma-optoAt(j)+1:ma-optoAt(j)+opto_lengths(j))=behExpts{j};
                withinExpt_cellAvs=cellAvs{j};
                allResponses(onCell:onCell+size(withinExpt_cellAvs,1)-1,ma-optoAt(j)+1:ma-optoAt(j)+opto_lengths(j))=withinExpt_cellAvs;
                onCell=onCell+size(withinExpt_cellAvs,1);
            end
        else
            allOptos=temp{1};
            allResponses=cellAvs{1};
            allBeh=tempbeh{1};
        end
        allCellResponses{i}=allResponses;
        allOptoData{i}=allOptos;
        allBehData{i}=allBeh;
        if i==1
            % sanity check plot allOptos
            figure();
            hold all;
            plot(allOptos');
        end
    end
    if any(diff(alltimes)~=0)
        error('sampling rate must be the same across expts');
    else
        returnTimes=0:alltimes(1):(size(allCellResponses{1},2)-1)*alltimes(1);
        timestep=alltimes(1);
    end
end






     