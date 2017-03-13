function data=fixRunningData(data)

threshForStep=0.05;
sameThresh=10^-3;

for i=1:size(data,1)
    currdata=data(i,:);
    flipInds=zeros(1,10);
%     disp(i);
    rawDiff=diff(currdata);
    flipInds=find(abs(rawDiff)>=threshForStep);
%     disp(length(flipInds));
    for j=1:length(flipInds)
        if flipInds(j)==1 
            if abs(currdata(flipInds(j)+1)-currdata(flipInds(j)+2))<sameThresh
                if flipInds(j)+1>length(currdata)
                else
                    currdata(flipInds(j)+1:flipInds(j)+2)=currdata(flipInds(j)+1:flipInds(j)+2)+currdata(flipInds(j))-currdata(flipInds(j)+1);
                end
            end
        elseif flipInds(j)+1>length(currdata) || flipInds(j)+2>length(currdata)
        elseif abs(currdata(flipInds(j)+1)-currdata(flipInds(j)+2))<sameThresh  
            currdata(flipInds(j)+1:flipInds(j)+2)=currdata(flipInds(j)+1:flipInds(j)+2)+currdata(flipInds(j))-currdata(flipInds(j)+1);
        elseif abs(currdata(flipInds(j))-currdata(flipInds(j)-1))<sameThresh
            currdata(flipInds(j)-1:flipInds(j))=currdata(flipInds(j)-1:flipInds(j))+currdata(flipInds(j)+1)-currdata(flipInds(j));
        end
    end
    data(i,:)=currdata;
end