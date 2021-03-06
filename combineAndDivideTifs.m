function combineAndDivideTifs(nSpatialPieces,nPixelsOverlap_x,nPixelsOverlap_y)
% Note that function assumes files sorted alphabetically by name is correct
% order
% nSpatialPieces is the number of pieces into which to divide the spatial
% extent of the tif movie, to the nearest matching option
% Options are 1,2,4,6,9,12 or 16 pieces
% nPixelsOverlap is number of pixels by which to overlap image pieces
% (e.g., if diameter of cell is about 12 pixels in Y dimension and 50
% pixels in X dimension, then set nPixelsOverlap_x=50 and
% nPixelsOverlap_y=12)

[file,pathname]=uigetfile('.tif','Choose which tif files to combine','MultiSelect','on');

try 
    % Movies can be different lengths
    combinedMovie=[];
    disp('Reading in tif files:');
    for i=1:length(file)
        disp(i);
        movie=tiffRead([pathname file{i}],'uint16');
        combinedMovie=cat(3,combinedMovie,movie);
    end
catch
    disp('Error in combining .tif files');
end

if ~isempty(combinedMovie)
    % Options
    segments=[1 2 4 6 9 12 16];
    if nSpatialPieces>16
        disp('Current max division is into 16 pieces');
        return
    elseif nSpatialPieces<1
        disp('nSpatialPieces must be greater than 0');
        return
    else
        temp=abs(segments-nSpatialPieces);
        [~,mi]=min(temp);
        nSpatialPieces=segments(mi);
    end
    if nSpatialPieces==1
        combined{1}=combinedMovie;
    elseif nSpatialPieces>1
        % Spatially sub-divide tif movie
        dimX=size(combinedMovie,2);
        dimY=size(combinedMovie,1);
        switch nSpatialPieces
            case 2
                combined=cell(1,2);
                rangeX1=1:floor(dimX/2 - nPixelsOverlap_x/2)+nPixelsOverlap_x;
                rangeX2=floor(dimX/2 - nPixelsOverlap_x/2):dimX;
                rangeY1=1:dimY;
                combined{1}=combinedMovie(rangeX1,rangeY1);
                combined{2}=combinedMovie(rangeX2,rangeY1);
            case 4
                combined=cell(1,4);
                rangeX1=1:floor(dimX/2 - nPixelsOverlap_x/2)+nPixelsOverlap_x;
                rangeX2=floor(dimX/2 - nPixelsOverlap_x/2):dimX;
                rangeY1=1:floor(dimY/2 - nPixelsOverlap_y/2)+nPixelsOverlap_y;
                rangeY2=floor(dimY/2 - nPixelsOverlap_y/2):dimY;
                combined{1}=combinedMovie(rangeY1,rangeX1,:);
                combined{2}=combinedMovie(rangeY1,rangeX2,:);
                combined{3}=combinedMovie(rangeY2,rangeX1,:);
                combined{4}=combinedMovie(rangeY2,rangeX2,:);
            case 6
                combined=cell(1,6);
                rangeX1=1:floor(dimX/3 - nPixelsOverlap_x/2)+nPixelsOverlap_x;
                rangeX2=floor(dimX/3 - nPixelsOverlap_x/2):floor(2*dimX/3 - nPixelsOverlap_x/2)+nPixelsOverlap_x;
                rangeX3=floor(2*dimX/3 - nPixelsOverlap_x/2):dimX;
                rangeY1=1:floor(dimY/2 - nPixelsOverlap_y/2)+nPixelsOverlap_y;
                rangeY2=floor(dimY/2 - nPixelsOverlap_y/2):dimY;
                combined{1}=combinedMovie(rangeY1,rangeX1,:);
                combined{2}=combinedMovie(rangeY1,rangeX2,:);
                combined{3}=combinedMovie(rangeY1,rangeX3,:);
                combined{4}=combinedMovie(rangeY2,rangeX1,:);
                combined{5}=combinedMovie(rangeY2,rangeX2,:);
                combined{6}=combinedMovie(rangeY2,rangeX3,:);
            case 9
                combined=cell(1,9);
                rangeX1=1:floor(dimX/3 - nPixelsOverlap_x/2)+nPixelsOverlap_x;
                rangeX2=floor(dimX/3 - nPixelsOverlap_x/2):floor(2*dimX/3 - nPixelsOverlap_x/2)+nPixelsOverlap_x;
                rangeX3=floor(2*dimX/3 - nPixelsOverlap_x/2):dimX;
                rangeY1=1:floor(dimY/3 - nPixelsOverlap_y/2)+nPixelsOverlap_y;
                rangeY2=floor(dimY/3 - nPixelsOverlap_y/2):floor(2*dimY/3 - nPixelsOverlap_y/2)+nPixelsOverlap_y;
                rangeY3=floor(2*dimY/3 - nPixelsOverlap_y/2):dimY;
                combined{1}=combinedMovie(rangeY1,rangeX1,:);
                combined{2}=combinedMovie(rangeY1,rangeX2,:);
                combined{3}=combinedMovie(rangeY1,rangeX3,:);
                combined{4}=combinedMovie(rangeY2,rangeX1,:);
                combined{5}=combinedMovie(rangeY2,rangeX2,:);
                combined{6}=combinedMovie(rangeY2,rangeX3,:);
                combined{7}=combinedMovie(rangeY3,rangeX1,:);
                combined{8}=combinedMovie(rangeY3,rangeX2,:);
                combined{9}=combinedMovie(rangeY3,rangeX3,:);
            case 12
                combined=cell(1,12);
                rangeX1=1:floor(dimX/4 - nPixelsOverlap_x/2)+nPixelsOverlap_x;
                rangeX2=floor(dimX/4 - nPixelsOverlap_x/2):floor(2*dimX/4 - nPixelsOverlap_x/2)+nPixelsOverlap_x;
                rangeX3=floor(2*dimX/4 - nPixelsOverlap_x/2):floor(3*dimX/4 - nPixelsOverlap_x/2)+nPixelsOverlap_x;
                rangeX4=floor(3*dimX/4 - nPixelsOverlap_x/2):dimX;
                rangeY1=1:floor(dimY/3 - nPixelsOverlap_y/2)+nPixelsOverlap_y;
                rangeY2=floor(dimY/3 - nPixelsOverlap_y/2):floor(2*dimY/3 - nPixelsOverlap_y/2)+nPixelsOverlap_y;
                rangeY3=floor(2*dimY/3 - nPixelsOverlap_y/2):dimY;
                combined{1}=combinedMovie(rangeY1,rangeX1,:);
                combined{2}=combinedMovie(rangeY1,rangeX2,:);
                combined{3}=combinedMovie(rangeY1,rangeX3,:);
                combined{4}=combinedMovie(rangeY1,rangeX4,:);
                combined{5}=combinedMovie(rangeY2,rangeX1,:);
                combined{6}=combinedMovie(rangeY2,rangeX2,:);
                combined{7}=combinedMovie(rangeY2,rangeX3,:);
                combined{8}=combinedMovie(rangeY2,rangeX4,:);
                combined{9}=combinedMovie(rangeY3,rangeX1,:);
                combined{10}=combinedMovie(rangeY3,rangeX2,:);
                combined{11}=combinedMovie(rangeY3,rangeX3,:);
                combined{12}=combinedMovie(rangeY3,rangeX4,:);
            case 16
                combined=cell(1,16);
                rangeX1=1:floor(dimX/4 - nPixelsOverlap_x/2)+nPixelsOverlap_x;
                rangeX2=floor(dimX/4 - nPixelsOverlap_x/2):floor(2*dimX/4 - nPixelsOverlap_x/2)+nPixelsOverlap_x;
                rangeX3=floor(2*dimX/4 - nPixelsOverlap_x/2):floor(3*dimX/4 - nPixelsOverlap_x/2)+nPixelsOverlap_x;
                rangeX4=floor(3*dimX/4 - nPixelsOverlap_x/2):dimX;
                rangeY1=1:floor(dimY/4 - nPixelsOverlap_y/2)+nPixelsOverlap_y;
                rangeY2=floor(dimY/4 - nPixelsOverlap_y/2):floor(2*dimY/4 - nPixelsOverlap_y/2)+nPixelsOverlap_y;
                rangeY3=floor(2*dimY/4 - nPixelsOverlap_y/2):floor(3*dimY/4 - nPixelsOverlap_y/2)+nPixelsOverlap_y;
                rangeY4=floor(3*dimY/4 - nPixelsOverlap_y/2):dimY;
                combined{1}=combinedMovie(rangeY1,rangeX1,:);
                combined{2}=combinedMovie(rangeY1,rangeX2,:);
                combined{3}=combinedMovie(rangeY1,rangeX3,:);
                combined{4}=combinedMovie(rangeY1,rangeX4,:);
                combined{5}=combinedMovie(rangeY2,rangeX1,:);
                combined{6}=combinedMovie(rangeY2,rangeX2,:);
                combined{7}=combinedMovie(rangeY2,rangeX3,:);
                combined{8}=combinedMovie(rangeY2,rangeX4,:);
                combined{9}=combinedMovie(rangeY3,rangeX1,:);
                combined{10}=combinedMovie(rangeY3,rangeX2,:);
                combined{11}=combinedMovie(rangeY3,rangeX3,:);
                combined{12}=combinedMovie(rangeY3,rangeX4,:);
                combined{13}=combinedMovie(rangeY4,rangeX1,:);
                combined{14}=combinedMovie(rangeY4,rangeX2,:);
                combined{15}=combinedMovie(rangeY4,rangeX3,:);
                combined{16}=combinedMovie(rangeY4,rangeX4,:);
            otherwise
                disp('Currently code supports only subdivision into 1,2,4,6,9,12 or 16 pieces');
        end
    end
    
    cropHere.rows=zeros(1,size(combined{1},1));
    cropHere.columns=zeros(1,size(combined{1},2));
    cropHere.columns(249:256)=1;
    cropHere.columns(1:18)=1;
    cropHere.rows(123:128)=1;
    combined{1}=cropThisMovie(combined{1},cropHere);
    disp('croppedMovie');
    combined{1}=uncropThisMovie(combined{1},cropHere);
    
    try
        for i=1:length(combined)
            tiffWrite(combined{i},['combined' num2str(i) '.tif'],pathname,'uint16');
        end
    catch
        % Sometimes, disk access fails due to intermittent
        % network problem. In that case, wait and re-try once:
        pause(60);
        for i=1:length(combined)
            tiffWrite(combined{i},['combined' num2str(i) '.tif'],pathname,'uint16');
        end
    end
else
    disp('Combined movie is empty');
end


function mov=cropThisMovie(mov,cropHere)

% Crop
mov=mov(:,cropHere.columns==0,:);
mov=mov(cropHere.rows==0,:,:);


function movStruct=uncropThisMovie(movStruct,cropHere)

persistent fakeData

% Uncrop, fill with nans
mov=nan(length(cropHere.rows),length(cropHere.columns),size(movStruct,3));
rowInds=find(cropHere.rows==0);
columnInds=find(cropHere.columns==0);
temp=nan(length(cropHere.rows),size(movStruct,2),size(movStruct,3));
for i=1:length(rowInds)
    ri=rowInds(i);
    temp(ri,:,:)=movStruct(i,:,:);
end
for i=1:length(columnInds)
    ci=columnInds(i);
    mov(:,ci,:)=temp(:,i,:);
end

if isempty(fakeData)
   fakeData=randi(2,1,size(mov,3)); 
else
   if length(fakeData)~=size(mov,3)
       fakeData=randi(2,1,size(mov,3)); 
   end
end

% Check that no pixels are always zero
% This check is required for movie to work with CNMF Ca2+ source detection
movsum=reshape(sum(mov,3),size(mov,1),size(mov,2));
[rowind,columnind]=find(movsum==0 | all(isnan(movsum),3));
mov(isnan(mov))=0;
for i=1:length(rowind)
    mov(rowind(i),columnind(i),:)=mov(rowind(i),columnind(i),:)+reshape(fakeData,1,1,length(fakeData));
end

movStruct=mov;