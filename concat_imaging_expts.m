function out=concat_imaging_expts(dirs)

% dirs is a list of directories containing expts' imaging data

for i=1:length(dirs)
    disp(i);
    currdir=dirs{i};
    listing=dir(currdir);
    if i==1
        % create new structure to contain all expts' data
        for j=1:length(listing)
            tempstruct=[];
            matind=regexp(listing(j).name,'.mat');
            if ~isempty(matind)
                % load all .mat files
                dataname=listing(j).name(1:matind-1);
                a=load([currdir '\' listing(j).name]);
                q=fieldnames(a);
                for k=1:length(q)
                    tempstruct.(q{k})=a.(q{k});
                end
                out.(dataname)=tempstruct.(dataname);
            end
        end
        out.cell_from_expt=cell(size(out.withinCellResponses));
        for j=1:length(out.cell_from_expt(1:end))
            out.cell_from_expt{j}=cell(size(out.withinCellResponses{j}));
            temp=out.cell_from_expt{j};
            for k=1:length(temp(1:end))
                temp{k}=i;
            end
            out.cell_from_expt{j}=temp;
        end
    else
        for j=1:length(listing)
            tempstruct=[];
            matind=regexp(listing(j).name,'.mat');
            if ~isempty(matind)
                % load all .mat files
                dataname=listing(j).name(1:matind-1);
                a=load([currdir '\' listing(j).name]);
                q=fieldnames(a);
                for k=1:length(q)
                    tempstruct.(q{k})=a.(q{k});
                end
                tempout.(dataname)=tempstruct.(dataname);
            end
        end
        tempout.cell_from_expt=cell(size(tempout.withinCellResponses));
        for j=1:length(tempout.cell_from_expt(1:end))
            tempout.cell_from_expt{j}=cell(size(tempout.withinCellResponses{j}));
            temp=tempout.cell_from_expt{j};
            for k=1:length(temp(1:end))
                temp{k}=i;
            end
            tempout.cell_from_expt{j}=temp;
        end
        % match experimental conditions
        if ~any(ismember(fieldnames(out),'condition_tags')) || ~any(ismember(fieldnames(tempout),'condition_tags'))
            error('condition_tags must be in these directories, so that can match experimental conditions across expts');
        end
        [list1,list2,list1_index_into_list2]=matchExptConditions(out.condition_tags,tempout.condition_tags);
        concatTheseConditions=list1_index_into_list2(~isnan(list1_index_into_list2));
        % match fields in tempout to fields in out
        fields_in_out=fieldnames(out);
        for j=1:length(fields_in_out)
            currf=fields_in_out{j};
            concatout.(currf)=[];
            if strcmp(currf,'times')
                currcat=[];
                part1=out.(currf);
                part2=tempout.(currf);
                if iscell(part1)
                    if isa(part2,'double')
                        currcat=part1;
                        currcat{length(part1)+1}=part2;
                    elseif iscell(part2)
                        currcat=[part1 part2];
                    else
                        error('unrecognized format for field times');
                    end
                elseif isa(part1,'double')
                    if isa(part2,'double')
                        currcat{1}=part1;
                        currcat{2}=part2;
                    elseif iscell(part2)
                        currcat{1}=part1;
                        currcat(2:2+length(part2)-1)=part2;
                    else
                        error('unrecognized format for field times');
                    end
                end     
            elseif strcmp(currf,'out')
                % discard
            elseif strcmp(currf,'trialByTrialBeh')
                % concatenate all trials, this will then match to condition
                % that includes all trials
                currcat=[];
                part1=out.(currf);
                part2=tempout.(currf);
                if iscell(part1)
                    if isa(part2,'double')
                        currcat=part1;
                        currcat{length(part1)+1}=part2;
                    elseif iscell(part2)
                        currcat=[part1 part2];
                    else
                        error('unrecognized format for field trialByTrialBeh');
                    end
                elseif isa(part1,'double')
                    if isa(part2,'double')
                        currcat{1}=part1;
                        currcat{2}=part2;
                    elseif iscell(part2)
                        currcat{1}=part1;
                        currcat(2:2+length(part2)-1)=part2;
                    else
                        error('unrecognized format for field trialByTrialBeh');
                    end
                end     
            elseif iscell(out.(currf))
                % make fields lists for match to list indices from matchExptConditions
                temp=out.(currf);
                out.(currf)=temp(1:end);
                temp=tempout.(currf);
                tempout.(currf)=temp(1:end);               
                concatInd=1;
                part1=out.(currf);
                part2=tempout.(currf);
                currcat=[];
                % concatenate experimental conditions that match, pulling
                % them to the top of structure
                if ~isempty(concatTheseConditions)
                    for k=1:length(list1_index_into_list2)
                        if isnan(list1_index_into_list2(k))
                        else
                            getThis1=part1{k};
                            getThis2=part2{list1_index_into_list2(k)};
                            if ~iscell(getThis1)
                                getThis1=part1(k);
                            end
                            if ~iscell(getThis2)
                                getThis2=part2(list1_index_into_list2(k));
                            end
                            if strcmp(currf,'condition_tags') % if these are the condition tags, take a consensus
                                currcat{concatInd}=getThis1; % will be the same tag
                            else % concatenate
                                currcat{concatInd}=[getThis1 getThis2];
                            end
                            concatInd=concatInd+1;
                        end
                    end
                end
                % add in experimental conditions from out that are not in
                % tempout
                notin=find(isnan(list1_index_into_list2));
                if ~isempty(notin)
                    for k=1:length(notin)
                        getThis1=part1{notin(k)};
                        if ~iscell(getThis1)
                            getThis1=part1(notin(k));
                        end
                        currcat{concatInd}=[getThis1];
                        concatInd=concatInd+1;
                    end
                end
                % add in experimental conditions from tempout that are not
                % in out
                condsintempout=1:length(part2);
                notin=condsintempout(~ismember(condsintempout,list1_index_into_list2));
                if ~isempty(notin)
                    for k=1:length(notin)
                        getThis2=part2{notin(k)};
                        if ~iscell(getThis2)
                            getThis2=part2(notin(k));
                        end
                        currcat{concatInd}=[getThis2];
                        concatInd=concatInd+1;
                    end
                end
            end
            concatout.(currf)=currcat;
        end
        out=concatout;
    end
end