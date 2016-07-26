function [Yr,b2,f2,Cn,Yk,Cf,Df,Ao]=readOrchestraOutput(fname)

fileId=fopen(fname);
if fileId==-1
    disp('Failed to open file');
end

gotYr=2;
gotb2=0;
gotf2=0;
gotCn=0;
gotYk=0;
gotCf=0;
gotDf=0;
gotAo=0;

safety=1;
Yr=[];
b2=[];
f2=[];
Cn=[];
Yk=[];
Cf=[];
Df=[];
Ao=[];
while any([gotYr~=2 gotb2~=2 gotf2~=2 gotCn~=2 gotYk~=2 gotCf~=2 gotDf~=2 gotAo~=2])
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
        b2=b2';
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
    elseif gotCn==2 && gotYk==0 && strcmp(tline,'    Yk')
        % Found start of Yk
        gotYk=1;
        disp('Found Yk');
    elseif gotCn==2 && gotYk==1 && ~isempty(tline)
        % Add to Yk
        Yk=[Yk; str2num(tline)];
    elseif gotCn==2 && gotYk==1 && isempty(tline)
        % End of Yk
        gotYk=2;
        disp('Finished Yk');
    elseif gotYk==2 && gotCf==0 && strcmp(tline,'    Cf')
        % Found start of Cf
        gotCf=1;
        disp('Found Cf');
    elseif gotYk==2 && gotCf==1 && ~isempty(tline)
        % Add to Cf
        Cf=[Cf; str2num(tline)];
    elseif gotYk==2 && gotCf==1 && isempty(tline)
        % End of Cf
        gotCf=2;
        disp('Finished Cf');
    elseif gotCf==2 && gotDf==0 && strcmp(tline,'    Df')
        % Found start of Df
        gotDf=1;
        disp('Found Df');
    elseif gotCf==2 && gotDf==1 && ~isempty(tline)
        % Add to Df
        Df=[Df; str2num(tline)];
    elseif gotCf==2 && gotDf==1 && isempty(tline)
        % End of Df
        gotDf=2;
        disp('Finished Df');        
    elseif gotDf==2 && gotAo==0 && strcmp(tline,'    Ao')
        % Found start of Ao
        gotAo=1;
        disp('Found Ao');
    elseif gotDf==2 && gotAo==1 && ~isempty(tline)
        % Add to Ao
        Ao=[Ao; str2num(tline)];
    elseif gotDf==2 && gotAo==1 && isempty(tline)
        % End of Ao
        gotAo=2;
        disp('Finished Ao');
        Ao=Ao';
    elseif gotYr==2 && gotb2==2 && gotf2==2 && gotCn==2 && gotYk==2 && gotCf==2  && gotDf==2  && gotAo==2
        % Got everything
        break
    end
    safety=safety+1;
    if safety>500000
        break
    end
end
        
    