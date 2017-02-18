function physData=findPhysData(obj,movieOrder,nameCommand)

% Look for phys command in same path as selected movie files
movName=obj.Movies{1};
dirBreaks=regexp(movName,'\','start');
physPath=movName(1:dirBreaks(end));
listing=dir([physPath nameCommand '*.mat']);
listingnames=cell(1,length(listing));
for i=1:length(listing)
    listingnames{i}=listing(i).name;
end

% Sort files by number
% movNumbers=nan(1,length(listingnames));
% for i=1:length(listingnames)
%     currMovName=listingnames{i};
%     currMovName=fliplr(currMovName);
%     startInd=regexp(currMovName,'\.','once');
%     currMovName=currMovName(startInd+1:end);
%     isNumberInd=regexp(currMovName,'\d');
%     numberMovName=currMovName(isNumberInd);
%     numberMovName=fliplr(numberMovName);
%     movNumbers(i)=str2num(numberMovName);
% end
% [~,inds]=sort(movNumbers);
% listingnames=listingnames(inds);

physData=[];
if ~isempty(listing)
    % Phys command saved during acquisition
    % For each movie file, get associated phys command
    
    for i=1:length(movieOrder)
        ci=movieOrder(i);
        s=obj.Movies{ci};
        parts=regexp(s,'\','split');
        movName=parts{end};
        fi=regexp(movName,'.tif');
        movName=movName(1:fi-1);
        movNumber=movName(end-2:end);
        % Check that expected file is in directory
        physFile=[nameCommand '_' num2str(str2num(movNumber)) '.mat'];
        if any(strcmp(listingnames,physFile))
            warning('off','MATLAB:unknownObjectNowStruct'); % import wave as struct
            % Load phys command
            w=load([physPath physFile]);
            warning('on','MATLAB:unknownObjectNowStruct'); 
            f=fieldnames(w);
            w=w.(f{1});
            if isempty(physData)
                physData=nan(length(movieOrder),length(w.data));
            end
            physData(i,:)=w.data;
        else
            disp(['Missing ' nameCommand ' file from movie directory']);
        end
    end
end

if isempty(physData)
    % No phys data 
    disp('No phys data found');
    times=[];
end
   
end