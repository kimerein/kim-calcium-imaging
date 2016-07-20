function [h,yOffsets]=plotCurvesOffset(h,data_x,data_y,sep,c,yOffsets)

if isempty(h)
    h=figure();
end

if isempty(yOffsets)
    yOffset=0;
    yOffsets=nan(length(data_x),1);
    for i=1:length(data_x)
        plot(data_x{i},data_y{i}+yOffset,'Color',c);
        hold on;
        yOffsets(i)=yOffset;
        yOffset=yOffset+max(data_y{i})+sep*max(data_y{i});
    end
else
    for i=1:length(data_x)
        plot(data_x{i},data_y{i}+yOffsets(i),'Color',c);
        hold on;
    end
end