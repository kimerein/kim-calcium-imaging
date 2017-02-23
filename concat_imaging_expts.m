function out=concat_imaging_expts(dirs)

% dirs is a list of directories containing expts' imaging data

for i=1:length(dirs)
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
                out.(dataname)=tempstruct;
            end
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
                tempout.(dataname)=tempstruct;
            end
        end
        disp('hi');
%         % match fields in tempout to fields in out
%         % will only continue w fields in both
%         fields_in_out=fieldnames(out);
%         for j=1:length(fields_in_out)
%             % concat
    end
end
            
            