function combineTifs()
% Note that function assumes files sorted alphabetically by name is correct
% order

[file,pathname]=uigetfile('.tif','Choose which tif files to combine','MultiSelect','on');

try
    % Assumes all movies are same size
    combinedMovie=[];
    for i=1:length(file)
        movie=tiffRead([pathname file{i}],'uint16');
        if i==1
            combinedMovie=nan([size(movie,1) size(movie,2) size(movie,3)*length(file)]);
        end
        combinedMovie(:,:,(i-1)*size(movie,3)+1:(i-1)*size(movie,3)+size(movie,3))=movie;
    end
catch
    try 
        % Movies different sizes
        combinedMovie=[];
        for i=1:length(file)
            movie=tiffRead([pathname file{i}],'uint16');
            combinedMovie=cat(3,combinedMovie,movie);
        end
    catch
        disp('Error in combining .tif files');
    end
end

if ~isempty(combinedMovie)
    try
        tiffWrite(combinedMovie,'combined.tif',pathname,'int16');
    catch
        % Sometimes, disk access fails due to intermittent
        % network problem. In that case, wait and re-try once:
        pause(60);
        tiffWrite(combinedMovie,'combined.tif',pathname,'int16');
    end
else
    disp('Combined movie is empty');
end