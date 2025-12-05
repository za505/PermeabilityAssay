%BacTrack.m
%Tracks bacterial growth from phase image stacks.  
%Customized for B. subtilis.

%INSTRUCTIONS FOR USE:
% Do NOT run all at once. Run ONE section at a time. 
%
%INPUT:
%basename: name of the image stack.
%dirname: the full pathname of the directory where you saved the (erased)
% image stack.
%lscale: microscope calibration in microns per pixels.
%sm: width of the Gaussian filter used in edge finder equals sm*sqrt(2).
%minL: minimum length of cells;
%minW: minimum width of cells;
%maxW: maximum width of cells;
%recrunch:0 or 1.  if you've already tracked the data set and just want to
%         re-plot the data enter 1.

%OUTPUT:
%T: number of time points.
%time: vector of length T with time points.
%ncells: number of individual cells tracked.
%lcell: ncells x T matrix of cell lengths.
%wcell: ncells x T matrix of cell widths.
%acell: ncells x T matrix of cell areas
%acell: ncells x T matrix of cell areas.
%B: ncells x T cell array with cell contours.
%mlines: ncells x T cell array with cell midlines

%Calls on the following m-files:
%norm16bit.m
%polefinder.m
%cellcurvature.m
%extrema.m
%effectivelength.m
%fig2pretty.m
%movingaverage.m

clear, close all

tic

%%%%%%%%%%%%%%%%%%%%%%%%%t%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%User Input
basename = '10012024_Exp2'; %Name of the image stack, used to save file.
dirname=[basename '/' basename '_phase/' basename '_erased/']; %Directory that the image stack is saved in.
savedir = ['~/tracking/']; %Directory to save the output .mat file to.

lscale=0.08; %Microns per pixel.
%tscale = [60, 300]; %in seconds, for 1 min --> 5 min frame rate movies
tscale = 60; % in seconds, for 1 min frame rate movies

thresh=0;%For default, enter zero.
IntThresh=10000;%Threshold used to enhance contrast. Default:35000
dr=1;%Radius of dilation before watershed default = 1
sm=2;%Parameter used in edge detection, default=2
minA=100;%Minimum cell area. default 50
maxA=2000; %maximum cell area. default 2000
cellLink=4;%Number of frames to ignore missing cells when tracking frame to frame
recrunch=0;%Display data from previously crunched data? 0=No, 1=Yes.
vis=0;%Display cell tracking? 0=No, 1=Yes.
checkhist=0;%Display image histogram? 0=No, 1=Yes.
troubleshooting=1; %pick the troubleshooting mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1
    cd(savedir)
    load([basename '_BT'])
    troubleshooting=0;
else

%Determine number of frames
cd(dirname);
curdir=cd;
directory=dir('*.tif');
T=length(directory);

% for uniform time scale
time = [0:T-1]*tscale;

% for two-time scale experiment
% tpoint1 = [0:11] * tscale(1);
% tpoint2 = tpoint1(end) + ([1:12] * tscale(2));
% time = [tpoint1, tpoint2];

cd(curdir);
path(dirname,path)

nc=zeros(1,T);
allcentroids=[];
cellnum=[];
tstamp=[];

%Pre-allocate matrices
acell=zeros(1,T);
wcell=zeros(1,T);
lcell=zeros(1,T);
DS=zeros(1,T);
B=cell(1,T);
pcell=zeros(1,T);
mlines=cell(1,T);
pxs = cell(1, T);

%Load first image
imagename=directory(1).name;
im=imread(imagename);
[imM,imN]=size(im);


%Track cells
for t=1:T
    t
    
    %Load image
    imagename=directory(t).name;

    im=imread(imagename);
    [imM,imN]=size(im);
    
    %De-speckle image
    im=medfilt2(im);
    
    %Normalize images
    ppix=0.5; %defauly 0.5
    im=norm16bit(im,ppix);
    
    %Enhance contrast
    imc=imcomplement(im);
    
    if checkhist==1 & (t >= T-5 | t <= 5);
        figure,imhist(imc),pause, close all
    end
    
    if thresh==0;
        [imcounts,bins]=imhist(imc);
        [imcounts,idx]=sort(imcounts);
        bins=bins(idx);
        thresh1=bins(end-1);
        if thresh1==65535
            thresh1=bins(end);
        end
    else
        thresh1=thresh;
    end

    imc=imadjust(imc,[thresh1/65535 1],[]);   

    %Find edges
    [ed2,thresh2]=edge(imc,'canny',[],sm*sqrt(2));

    %Clean image
    cc=bwconncomp(ed2,8); %default 8
    stats=regionprops(cc,imc,'Area','MeanIntensity');
    idx=find([stats.Area]>minA&[stats.Area]<maxA&[stats.MeanIntensity]>IntThresh);
    ed2=ismember(labelmatrix(cc),idx);
        
    %Close gaps in edges
    despurred=bwmorph(ed2,'spur');
    spurs=ed2-despurred;
    [spy,spx]=find(spurs);
    for k=1:length(spx)
        ed2(spy(k)-1:spy(k)+1,spx(k)-1:spx(k)+1)=ed2(spy(k)-1:spy(k)+1,spx(k)-1:spx(k)+1)+rot90(ed2(spy(k)-1:spy(k)+1,spx(k)-1:spx(k)+1),2);
        ed2(spy(k),spx(k))=1;
    end
    ed2=bwmorph(ed2,'bridge'); 
            
    se=strel('disk',dr);
    ed2=imdilate(ed2,se);
    ed2=bwmorph(ed2,'thin',2);
    
    %Identify cells based on size and intensity
    ed3=~ed2;
    ed3(1,:)=ones(1,imN);
    ed3(end,:)=ones(1,imN);
    ed3(:,1)=ones(imM,1);
    ed3(:,end)=ones(imM,1);
              
    cc=bwconncomp(ed3,4);
    stats=regionprops(cc,imc,'Area','MeanIntensity');

    % idx=find([stats.Area]>minA&[stats.Area]<maxA&[stats.MeanIntensity]>3e4);
    idx=find([stats.Area]>minA&[stats.Area]<maxA&[stats.MeanIntensity]>1e4);
    ed4=ismember(labelmatrix(cc),idx);

    %Find cell areas and centroids
    bw=bwmorph(ed4,'thicken');
    [P,bw]=bwboundaries(bw,4,'noholes');
    stats=regionprops(bw, 'Centroid','PixelIdxList');
    
    L=bwlabel(bw);    
    %labels(:,:,t)=L;
    %labels2(:,:,t)=bw;

    nc(t)=length(P);
    cents=cat(1,stats.Centroid);
    centroids=cents;
     
    allcentroids=[allcentroids;centroids];
    tstamp=[tstamp;ones(nc(t),1)*t];
    cellnum=[cellnum;(1:nc(t))'];

    %save pixel id list
    for n=1:nc(t)
        pxls{n,t}=stats(n).PixelIdxList;
    end

    pxs{t} = pxls(:, t);

    if vis==1 & (t >= T-5 | t <= 5)
       figure
       imshow(bw)
       hold on

       scatter(centroids(:,1), centroids(:,2),'*r')

      pause
      close all
    end
    
    %toc

end

end

["Finished Part 1"]

%%
%Track cells frame to frame
tracks=zeros(size(im));
rcents=round(allcentroids);
linind=sub2ind(size(im),rcents(:,2),rcents(:,1));
tracks(linind)=1;

nhood=[0,1,0;1,1,1;0,1,0];
tracks=imdilate(tracks,strel('disk',cellLink));
overlay1=imoverlay(im,tracks,[.3 1 .3]);

[tracksL,ncells]=bwlabel(tracks,8);

pixs=cell(ncells,T);
lcents=length(allcentroids);

for i=1:lcents
    cellid=tracksL(linind(i));
    pixs{cellid,tstamp(i)}=pxs{tstamp(i)}{cellnum(i), 1};
end

for i=1:T
    for j=1:ncells
        if isempty(pixs{j,i})==0 
            cell_id_na(j,i)=j; % cell id not aligned/ not lineage tracked
        end
    end
end

for i = flip(1:T)
    for j=1:ncells
        if isempty(pixs{j,i})==1
            pixs{j,i} = [NaN];
        end
    end
end

%Align cells from frame to frame, correct pixs to track lineage
for i=T 
    for j=1:ncells
        cell_id(j,i)=j; %makes cell id vector
        pixels{j,i}=pixs{j,i};%Seed pixels
    end
end

%Tracks cells backward through time (note: if a cell isn't tracked in the final frame, it's not tracked for the entire movie) 
for i=flip(2:T)
    for j=1:ncells
        if ~isnan(pixels{j,i}) %if nonNaN
            nidx=cellfun('length', pixels(j,:));
            nnidx=nidx>1;
            ncols=min(find(nnidx)); %last nonNaN value
            Lia1=ismember(pixs{j,i-1},pixels{j,ncols},'rows'); %look within that array's row
            if sum(Lia1)>1
                cell_id(j,i-1)=j;
                pixels{j,i-1}=pixs{j,i-1};
            else
                oidx = find(cellfun(@(qset) sum(ismember(pixels{j,i}, qset)), pixs(:, i-1))>1); %look across all rows within the array
                if length(oidx)>1
                    cell_id(j,i-1)=j;
                    pixels{j, i-1} = unique(cell2mat(pixs(oidx, i-1)));
                else
                    cell_id(j,i-1)=NaN;
                    pixels{j,i-1}=NaN;
                end
            end
        else %if the value is NaN
            index=cellfun('length', pixels(j,:));
            nindex=index>1;
            ncols=min(find(nindex));
            Lia1=ismember(pixels{j,ncols}, pixs{j,i-1},'rows'); %find overlap within the array's row
            if sum(Lia1)>1
                cell_id(j,i-1)=j;
                pixels{j,i-1}=pixs{j,i-1};
            else
                cell_id(j,i-1)=NaN;
                pixels{j,i-1}=NaN;
            end
        end
    end
end

%how many time points are untracked?
lga = cellfun(@(A) sum(isnan(A)), pixels);
sia = sum(lga, 2);
eidx = find(sia >= T-2); %include cells with at least three time points

% Fill NaNs tracking backward
for i = flip(1:T)
    for j=setdiff(1:ncells, eidx)
        if isnan(pixels{j,i}) %if a value is NaN
            index=cellfun('length', pixels(j, 1:T-2));
            nindex=index==1;
            ncols=max(find(nindex)); %look for the last NaN value
            if isempty(ncols)
                ncols = i;
            end
            oidx = find(cellfun(@(qset) sum(ismember(pixels{j,ncols+1}, qset)), pixs(:, i))>1); %find pixel overlaps with the first nonNaN value after the last NaN value
            if length(oidx)>1
                cell_id(j,i)=j;
                pixels{j, i} = unique(cell2mat(pixs(oidx, i)));
            elseif length(oidx)==1
                cell_id(j,i)=j;
                pixels{j, i} = pixs{oidx, i};
            else
                cell_id(j,i)=NaN;
                pixels{j,i}=NaN;
            end
        end
    end
end

%how many time points are untracked? remove these cells
lga = cellfun(@(A) sum(isnan(A)), pixels);
sia = sum(lga, 2);
fidx = find(sia < T-2); %include cells with at least three time points

pixels = pixels(fidx, :);

[ncells, ~] = size(pixels);
pxls = cell(ncells, T);
B = cell(ncells, T);

%label = zeros(imM, imN, T);
%labels = zeros(imM, imN, T);

for t=1:T

    imBW = zeros(size(im));
    imL = zeros(size(im));

    for n=1:ncells
        if ~isnan(pixels{n,t})
            imB = zeros(size(im));
            imB(pixels{n,t})=1;

            imBW(pixels{n,t})=1;
            imL(pixels{n,t})=n;

            [P,bw]=bwboundaries(imB,4,'noholes');
            stats=regionprops(bw,'Area','Centroid','PixelIdxList');

            if length(P)==1

                acell(n,t)=[stats.Area]';

                rP=[P{1}(:,2),P{1}(:,1)];
                px=[rP(1:end-1,1);rP(1:end-1,1);rP(:,1)];
                py=[rP(1:end-1,2);rP(1:end-1,2);rP(:,2)];
                sp=length(rP);
                dS=sqrt(diff(px).^2+diff(py).^2);
                S=[0 cumsum(dS)'];

                px=csaps(S,px,0.05,S);
                py=csaps(S,py,0.05,S);

                px=px(sp+1:2*sp);
                py=py(sp+1:2*sp);

                px(end)=px(1);
                py(end)=py(1);

                dS=sqrt(diff(px).^2+diff(py).^2);
                S=[0 cumsum(dS)];
                ls=length(S);
                DS(n,t)=S(end)/(ls-1);
                Sn=(0:DS(n,t):S(end));
                nx=spline(S,px,Sn);
                ny=spline(S,py,Sn);

                B{n,t}=[nx',ny'];
                pxls{n,t}=stats.PixelIdxList;

                X=B{n,t}(:,1);
                Y=B{n,t}(:,2);   

                [sX,~]=size(X);

                %Find poles
                [X,Y,pcell(n,t)]=polefinder(X,Y);

                %Create mesh
                npts=min(pcell(n,t),sX-pcell(n,t)+1);
                S=(0:DS(n,t):(sX-1)*DS(n,t));

                s1=(0:S(pcell(n,t))/(npts-1):S(pcell(n,t)));
                s2=(S(pcell(n,t)):(S(end)-S(pcell(n,t)))/(npts-1):S(end));
                xc1=spline(S(1:pcell(n,t)),X(1:pcell(n,t)),s1);
                xc2=spline(S(pcell(n,t):end),X(pcell(n,t):end),s2);
                yc1=spline(S(1:pcell(n,t)),Y(1:pcell(n,t)),s1);
                yc2=spline(S(pcell(n,t):end),Y(pcell(n,t):end),s2);
                xc2=fliplr(xc2);
                yc2=fliplr(yc2);

                %Calculate midline
                mlines{n,t}=[(xc1+xc2)'/2,(yc1+yc2)'/2];
                dxy=diff(mlines{n,t}).^2;
                dl=sqrt(dxy(:,1)+dxy(:,2));
                lcell(n,t)=sum(dl);

                %Calculate width
                ls=[0 cumsum(dl)'];
                [~,mpos1]=min(abs(ls/lcell(n,t)-0.25));
                [~,mpos2]=min(abs(ls/lcell(n,t)-0.75));

                widths=sqrt((xc1-xc2).^2+(yc1-yc2).^2);
                wcell(n,t)=(widths(mpos1)+widths(mpos2))/2;
        end
      end
    end

    %label(:, :, t) = imL;
    %labels(:, :, t) = imBW;

end

lcell(lcell==0)=NaN;
wcell(wcell==0)=NaN;
acell(acell==0)=NaN;
pcell(pcell==0)=NaN;

%Dimsionalize the variables
lcell=lcell*lscale;
wcell=wcell*lscale;
acell=acell*lscale^2;

%Calculate total length of all cells
wcell(isnan(wcell))=0; 
lcell(isnan(lcell))=0;
acell(isnan(acell))=0;
for t=1:T
    ltotal(t)=sum(nonzeros(lcell(:,t)));
    atotal(t)=sum(nonzeros(acell(:,t)));
end

% Set zero values back to NaN
lcell(lcell==0)=NaN;
wcell(wcell==0)=NaN;
acell(acell==0)=NaN;

["Finished Tracking"]


%% Troubleshooting, Manually checking the tracking with the final image
if troubleshooting==1 % check cell by cell across all time points
    for k=1:ncells
       k

       widthInches = 8;
       heightInches = 6;
       figure('Units', 'inches', 'Position', [1, 1, widthInches, heightInches])
       imshow(im, [])
       hold on

       for t=1:T
         if isempty(B{k,t})==0
            plot(B{k,t}(:,1),B{k,t}(:,2),'-r')
         else
             continue
         end
       end

      pause
      close all
    end

elseif troubleshooting==2 % check all cells at each time point

    cd(dirname)

   for t=1:T
   t

    %Load image
    imagename=directory(t).name;
    im=imread(imagename);

   figure
   imshow(im)
   hold on

   for k=1:ncells
     if isempty(B{k,t})==0
        plot(B{k,t}(:,1),B{k,t}(:,2),'-r')
     else
         continue
     end
   end

  pause
  close all

   end 

elseif troubleshooting==3 % check cell by cell across select time points

     for k = 1:ncells
            k
            
            %Load image
            imagename=directory(t).name;
            im=imread(imagename);


            figure
            imshow(im, [])
            hold on

           for t=[1, 5, 25, 30, 40, 10, 30, 40, 50]
                if isempty(B{k,t})==0
                    plot(B{k,t}(:,1),B{k,t}(:,2),'-g')              
                end
           end
          pause
          close all
     end
end

%% Manually check cells
set = []; % which cells do you want to see?
plot_with_visibility_panel(time./60,lcell(set,:));
xlabel('Time (min)')
ylabel('Length (\mum)')
fig2pretty
xline(4, '--k')

%% Manually filter cells
%Filter out cells found at only one or two time points
kcells = [1, 2, 3, 5, 7, 9, 16]; %cells to keep

idx = ismember(1:ncells, kcells);
lcell = lcell(idx,:);
wcell=wcell(idx,:);
acell = acell(idx,:);
pcell = pcell(idx,:);
B = B(idx,:);
pixels = pixels(idx,:);
mlines = mlines(idx,:);
[ncells,~]=size(lcell);

cd(savedir)
save([basename '_BT'])


%% save plot

%cd(savedir)
% figure(1), title('Cell Length vs. Time'), hold on
% for i=1:ncells  
%     plot(time./60,lcell(i,:))
% end
% xlabel('Time (min)')
% ylabel('Length (\mum)')
% fig2pretty
% xline(5, '--k')
% saveas(gcf,[basename,'_lTraces.png'])
% saveas(gcf,[basename,'_lTraces.fig'])

%% Function
function plot_with_visibility_panel(x, data)
    % Create figure
    f = figure('Name', 'Plot with Visibility Panel', 'Position', [100 100 800 600]);

    % Plot each row and store handles
    nrows = size(data, 1);
    ax = axes('Parent', f, 'Position', [0.35, 0.1, 0.6, 0.85]);
    hold(ax, 'on');
    lineHandles = gobjects(nrows, 1);

    for i = 1:nrows
        lineHandles(i) = plot(ax, x, data(i, :), 'DisplayName', sprintf('Row %d', i));
    end

    legend(ax, 'show');

    % Create panel with checkboxes
    panel = uipanel('Parent', f, 'Title', 'Toggle Visibility', ...
        'Position', [0.02 0.1 0.3 0.85]);

    for i = 1:nrows
        uicontrol('Style', 'checkbox', ...
            'Parent', panel, ...
            'String', sprintf('Row %d', i), ...
            'Value', 1, ...
            'Position', [10, 400 - 25*i, 100, 20], ...
            'Callback', @(src, ~) toggle_visibility(src, lineHandles(i)));
    end
end

function toggle_visibility(checkbox, line)
    if checkbox.Value
        line.Visible = 'on';
    else
        line.Visible = 'off';
    end
end


