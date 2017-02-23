function [list1,list2,list1_index_into_list2]=matchExptConditions(tags1,tags2)

% makes a reference structure indicating how to concatenate data from more
% than one expt
list1=tags1(1:end);
list2=tags2(1:end);

% find how these match up
list1_index_into_list2=nan(1,length(list1));
for i=1:length(list1)
    currstr=list1{i};
    for j=1:length(list2)
        if strcmp(currstr,list2{j})==1
            list1_index_into_list2(i)=j;
        end
    end
end

% Make sure that lists are column vectors
list1=reshape(list1,length(list1),1);
list2=reshape(list2,length(list2),1);