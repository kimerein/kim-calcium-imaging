function [allResponses,optoTiming]=align_cellResponses(alldata,doAverage)

if doAverage==1
    % cycle through expt conditions
    for i=1:length(alldata.optoForProfile)
        temp=alldata.optoForProfile{i};
        optoExpts=cell(1,length(temp));
        exptTimes=cell(1,length(temp));
        cellAvs=cell(1,length(temp));
        for j=1:length(temp)
            optoStims=temp(j);
            av=nanmean(optoStims,1);
            optoExpts{j}=av;
            
        
    
    temp=testout.optoForProfile{1};
    temp1=temp(1);
    temp11=temp1{1};
    figure(); plot(nanmean(temp11,1));