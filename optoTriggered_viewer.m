function optoTriggered_viewer(C_df,Cn,f2,Ao,b2,options,acrossTrials,opto,beh,times)

% Get options used for Ca2+ source extraction
defoptions = CNMFSetParms;
if nargin < 7 || isempty(options); options = []; end
if ~isfield(options,'d1') || isempty(options.d1); d1 = input('What is the total number of rows? \n'); else d1 = options.d1; end          % # of rows
if ~isfield(options,'d2') || isempty(options.d2); d2 = input('What is the total number of columns? \n'); else d2 = options.d2; end          % # of columns
% if ~isfield(options,'normalize') || isempty(options.normalize); options.normalize = ones(size(Ao,1),1); end
%     sn = options.normalize;
if ~isfield(options,'plot_df') || isempty(options.plot_df); options.df = defoptions.plot_df; end
plot_df = options.plot_df;
if ~isfield(options,'make_gif') || isempty(options.make_gif); options.make_gif = defoptions.make_gif; end
make_gif = options.make_gif;
if ~isfield(options,'save_avi') || isempty(options.save_avi); options.save_avi = defoptions.save_avi; end
save_avi = options.save_avi;
if ~isfield(options,'sx') || isempty(options.sx); options.sx = defoptions.sx; end
sx = min([options.sx,floor(d1/2),floor(d2/2)]);
%if ~isfield(options,'pause_time') || isempty(options.pause_time); options.pause_time = defoptions.pause_time; end
%pause_time = options.pause_time;
if isfield(options,'name') && ~isempty(options.name);
    name = [options.name,'_components'];
else
    name = [defoptions.name,'_components'];
end

% Count
nr = size(C_df,1);     % number of ROIs
nb = size(f2,1);     % number of background components

% Video save settings
if save_avi
    vidObj = VideoWriter([name,'.avi']);
    set(vidObj,'FrameRate',1);
    open(vidObj);
end
thr = 0.95;
fig = figure('Visible','off');
set(gcf,'Position',2*[300,300,960,480]);
set(gcf,'PaperPosition',2*[300,300,960,480]);
int_x = zeros(nr,2*sx);
int_y = zeros(nr,2*sx);
cm = com(Ao,d1,d2);

% Create a figure and axes
% ax = axes('Units','DF/F');

% Create slider
sld = uicontrol('Style', 'slider',...
    'Min',1,'Max',nr+nb,'Value',1,'SliderStep',[1/(nr+nb-1) 1],...
    'Position', [150 20 800 20],...
    'Callback', @surfzlim);

% Add a text uicontrol to label the slider.
txt = uicontrol('Style','text',...
    'Position',[400 45 120 20],...
    'String','Component');

% Make figure visble after adding all components
fig.Visible = 'on';
plot_component(1)

% This code uses dot notation to set properties.
% Dot notation runs in R2014b and later.
% For R2014a and earlier: set(f,'Visible','on');


    function surfzlim(source,callbackdata)
        i = source.Value;
        plot_component(round(i))
        % For R2014a and earlier:
        % i = get(source,'Value');
        
        
%         if save_avi
%             currFrame = getframe(fig);
%             writeVideo(vidObj,currFrame);
%         else
%             pause(0.05);
%         end
    end

    function plot_component(i)
       if i <= nr
            subplot(3,2,5);
            Atemp = reshape(Ao(:,i),d1,d2);
            int_x(i,:) = round(cm(i,1)) + (-(sx-1):sx);
            if int_x(i,1)<1
                int_x(i,:) = int_x(i,:) + 1 - int_x(i,1);
            end
            if int_x(i,end)>d1
                int_x(i,:) = int_x(i,:) - (int_x(i,end)-d1);
            end
            int_y(i,:) = round(cm(i,2)) + (-(sx-1):sx);
            if int_y(i,1)<1
                int_y(i,:) = int_y(i,:) + 1 - int_y(i,1);
            end
            if int_y(i,end)>d2
                int_y(i,:) = int_y(i,:) - (int_y(i,end)-d2);
            end
            Atemp = Atemp(int_x(i,:),int_y(i,:));
            imagesc(int_x(i,:),int_y(i,:),Atemp); axis square;
        end
        subplot(3,2,[1,3]);
        if i <= nr
            cla
            imagesc(2*Cn); axis equal; axis tight; axis off; hold on;
            A_temp = full(reshape(Ao(:,i),d1,d2));
            A_temp = medfilt2(A_temp,[3,3]);
            A_temp = A_temp(:);
            [temp,ind] = sort(A_temp(:).^2,'ascend');
            temp =  cumsum(temp);
            ff = find(temp > (1-thr)*temp(end),1,'first');
            if ~isempty(ff)
                [~,ww] = contour(reshape(A_temp,d1,d2),[0,0]+A_temp(ind(ff)),'LineColor','k');
                ww.LineWidth = 2;
            end
            title(sprintf('Component %i ',i),'fontsize',16,'fontweight','bold'); drawnow; %pause;
        else
            cla
            imagesc(reshape(b2(:,i-nr),d1,d2)); axis equal; axis tight;
            title('Background component','fontsize',16,'fontweight','bold'); drawnow;
        end
        subplot(3,2,[2,4,6]);
        if i <= nr
%             plot(1:T,Y_r(i,:)/Df(i),'linewidth',2); hold all; 
            currAv=acrossTrials{i};
            plot(times,nanmean(currAv,1),'linewidth',2,'Color','r');
            hold all;
%             Overlay opto stim times
            plot(times,(opto./max(opto))*max(nanmean(currAv,1)),'linewidth',2);
            plot(times,beh,'linewidth',1);
            if plot_df
                title(sprintf('Component %i (calcium DF/F value)',i),'fontsize',16,'fontweight','bold');
            else
                title(sprintf('Component %i (calcium raw value)',i),'fontsize',16,'fontweight','bold');
            end
            leg = legend('Inferred','Opto Stim','Angular Velocity');
            set(leg,'FontSize',14,'FontWeight','bold'); 
            plot(times,nanmean(currAv,1)+nanstd(currAv,[],1)./sqrt(size(currAv,1)),'linewidth',1,'Color','r');
            plot(times,nanmean(currAv,1)-nanstd(currAv,[],1)./sqrt(size(currAv,1)),'linewidth',1,'Color','r');
            xlabel('Time (s)');
            drawnow;
            hold off;
            if make_gif
%                 frame = getframe(fig); %getframe(1);
%                 im = frame2im(frame);
%                 [imind,clm] = rgb2ind(im,256);
%                 if i == 1;
%                     imwrite(imind,clm,[name,'.gif'],'gif', 'Loopcount',inf);
%                 else
%                     imwrite(imind,clm,[name,'.gif'],'gif','WriteMode','append');
%                 end
            else
%                 if i < nr+nb && ~save_avi
%                     fprintf('component %i. Press any key to continue.. \n', i);
%                     if pause_time == Inf;
%                         pause;
%                     else
%                         pause(pause_time);
%                     end

%                 end
            end
        else
            plot(1:T,f2(i-nr,:)); title('Background activity','fontsize',16,'fontweight','bold');
            drawnow;
            if make_gif
                frame = getframe(fig); %getframe(1);
                im = frame2im(frame);
                [imind,clm] = rgb2ind(im,256);
                if i == 1;
                    imwrite(imind,clm,[name,'.gif'],'gif', 'Loopcount',inf);
                else
                    imwrite(imind,clm,[name,'.gif'],'gif','WriteMode','append');
                end
            else
%                 if i < nr+nb && ~save_avi
%                     fprintf('background component %i. Press any key to continue.. \n', i-nr);
% %                     if pause_time == Inf;
% %                         pause;
% %                     else
% %                         pause(pause_time);%                     end
%                 end
            end
        end 
    end
end