function cell_by_traces=convert_trials_to_continuous(traces)

for i=1:size(traces,1)
    temp=[];
    for j=1:size(traces,2)
        temp=[temp traces{i,j}];
    end
    if i==1
        cells_by_traces=nan(size(traces,1),length(temp));
        cell_by_traces(i,:)=temp;
    else
        cell_by_traces(i,:)=temp;
    end
end

    