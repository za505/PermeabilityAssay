%DiffusionMeasure.m
%Author: Zarina Akbary
%Calculates fluor. intensity over time. Incorporates BTfluo.m code.

clear, close all

%INSTRUCTIONS FOR USE:
%run BacTrack.m first

%INPUT
%basename: experiments of interest
%dirname: where .mat files are stored
%channels: list of directories containing fluorescent image stacks to quantify.

%OUTPUT:
%icell_intensity = cell x time matrix of fluoresc. intensities
%bg_intensity = 1 x time matrix of mean background intensities

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename = '10012024_Exp2'; %Name of the image stack, used to save file.

dirname= ['~/tracking/']; %Directory that the BT.mat files is saved in
savedir = ['~/fluorescence'];  %Directory to save the output .mat file to.
channels={[basename '/' basename '_mNeonGreen/' basename '_aligned/']}; %Directory that the image stack is saved in.

recrunch=0;
troubleshoot=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if recrunch==1

    cd(savedir)
    load([basename '_dm2.mat'])
    
    % omit outliers/mistracked cells that were missed in BacTrack.m
    delind = [13];

    idx=setdiff(1:ncells, delind);
    icell_intensity=icell_intensity(idx,:);
    sum_intensity=sum_intensity(idx,:);
    lcell=lcell(idx,:);
    wcell=wcell(idx,:);
    acell=acell(idx,:);
    vlme=vlme(idx,:);
    B=B(idx,:);
    ncells=length(idx);
    pixels=pixels(idx,:);
    pixel_num = pixel_num(idx,:);
    ncells = length(idx);
    
else 
    
    for i=1:length(channels)
        cd(channels{i}); 
        fluo_directory{i}=dir('*.tif');
    end
    
    %go to directory where .mat files are stored
    cd(dirname)
    
    load([basename '_BT'], 'B', 'T', 'ncells', 'time', 'pixels', 'lcell', 'wcell', 'acell')
       
    bg_intensity=nan(1,T);
    time=time./60;
   
    r = wcell./2;
    h = lcell-wcell;
    vlme = (pi * r.^2 .* h) + (4/3 * pi .* r.^3);

    for i=1:length(channels)
        
        cd(channels{i});
        
        % select the background region (pre-lysis, exclude region with
        % cells)
        imagename=fluo_directory{i}(5).name;
        [p1, p2]=getBackground(imagename);
        close all
     
        for t=1:T
            t
            imagename=fluo_directory{i}(t).name;
            im=imread(imagename);
            
            % calculate the mean fluorescence of the background region
            bg_intensity(1, t)=measureBackground(im, p1, p2);   

            % calculate the mean pixel intensity within the cell
           for n=1:ncells
            if ~isnan(pixels{n,t})
                icell_intensity(n,t)=mean(im(pixels{n,t})); 
                sum_intensity(n,t)= sum(im(pixels{n,t}));
                pixel_num(n, t) = length(pixels{n,t});
            else
                icell_intensity(n,t)=NaN;
                sum_intensity(n,t)= NaN;
                pixel_num(n, t) = NaN;
            end
           end

        end  
       
    end
end


% Plot data
cd(savedir)

%plot cellular fluorescence traces
figure, hold on
for i = 1:ncells
    plot(time, icell_intensity(i, :), '-g')
end
plot(time, bg_intensity, '-b')
title('Intensity vs Time')
%xline(4, '--k') % start of detergent perfusion
%xline(7, '--k') % end of detergent perfusion
xlabel('Time (min)')
ylabel('Cellular Intensity (A.U.)')
ylim([0 Inf])
%saveas(gcf, [basename,'_cellIntensity_dm.fig'])
%saveas(gcf, [basename,'_cellIntensity_dm.png'])

%plot cellular fluorescence traces
figure, hold on
for i = 1:ncells
    ycell = icell_intensity(i, :);
    y = ycell - bg_intensity;
    x = time;
    plot(x, y, '-r')
end
title('Cell Intensity vs Time')
xlabel('Time (min)')
ylabel('Cellular Intensity - Background Intensity (A.U.)')


cd(savedir)
save([basename '_dm.mat'])


%% Functions
 function [p1, p2]=getBackground(imagename)

         %Load last image
         im2=imread(imagename);

         %Determine Background
         figure,imshow(im2,[]), hold on, title('Select Background')
         %figure,imshow(im2,[50 5000]), hold on, title('Select Background')
         k=waitforbuttonpress;
         set(gcf,'Pointer')
         hold on
         axis manual
         point1=get(gca,'CurrentPoint');
         finalRect=rbbox;
         point2=get(gca,'CurrentPoint');
         point1=point1(1,1:2);
         point2=point2(1,1:2);
         point1(point1<1)=1;
         point2(point2<1)=1;
         p1=min(point1,point2);%Calculate locations
         p2=max(point1,point2);
         offset = abs(point1-point2);%And dimensions
         x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
         y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
         plot(x,y);
         p1=round(p1);
         p2=round(p2);  
 end 

 function bglevel = measureBackground(im2, p1, p2)

         %Determine background
         backim=im2(p1(2):p2(2),p1(1):p2(1));
         bglevel=mean(mean(backim));

 end 


