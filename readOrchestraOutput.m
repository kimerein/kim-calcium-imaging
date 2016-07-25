function [Yr, b2, f2, Cn]=readOrchestraOutput(fname)

fileId=fopen(fname);
if fileId==-1
    disp('Failed to open file');
end

gotYr=0;
gotb2=0;
gotf2=0;
gotCn=0;

safety=1;
Yr=[];
b2=[];
f2=[];
Cn=[];
while any([gotYr==0 gotb2==0 gotf2==0 gotCn==0])
    if mod(safety,1000)==1
        disp(safety);
    end
    tline=fgetl(fileId);
    if gotYr==0 && strcmp(tline,'    Yr')
        % Found start of Yr
        gotYr=1;
        disp('Found Yr');
    elseif gotYr==1 && ~isempty(tline)
        % Add to Yr
        Yr=[Yr; str2num(tline)];
    elseif gotYr==1 && isempty(tline)
        % End Yr
        gotYr=2;
        disp('Finished Yr');
    elseif gotYr==2 && gotb2==0 && strcmp(tline,'    b2')
        % Found start of b2
        gotb2=1;
        disp('Found b2');
    elseif gotYr==2 && gotb2==1 && ~isempty(tline)
        % Add to b2
        b2=[b2; str2num(tline)];
    elseif gotYr==2 && gotb2==1 && isempty(tline)
        % End of b2
        gotb2=2;
        disp('Finished b2');
    elseif gotb2==2 && gotf2==0 && strcmp(tline,'    f2')
        % Found start of f2
        gotf2=1;
        disp('Found f2');
    elseif gotb2==2 && gotf2==1 && ~isempty(tline)
        % Add to f2
        f2=[f2; str2num(tline)];
    elseif gotb2==2 && gotf2==1 && isempty(tline)
        % End of f2
        gotf2=2;
        disp('Finished f2');
    elseif gotf2==2 && gotCn==0 && strcmp(tline,'    Cn')
        % Found start of Cn
        gotCn=1;
        disp('Found Cn');
    elseif gotf2==2 && gotCn==1 && ~isempty(tline)
        % Add to Cn
        Cn=[Cn; str2num(tline)];
    elseif gotf2==2 && gotCn==1 && isempty(tline)
        % End of Cn
        gotCn=2;
        disp('Finished Cn');
    elseif gotYr==2 && gotb2==2 && gotf2==2 && gotCn==2
        % Got everything
        break
    end
    safety=safety+1;
    if safety>500000
        break
    end
end
        
    