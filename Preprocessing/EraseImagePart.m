%EraseImagePart.m
%Rico Rojas, edited by Zarina Akbary
%This routine lets the user select portions of an image stack to erase 
%(make background).   

%INSTRUCTIONS FOR USE:
%Duplicate the directory with the aligned phase images. When you run the
%program, the image specified in pickRegion will open. Select the regions
%you want to delete with the cursor and then press Enter. The program
%writes over the original image stack, so if you want a backup stack, save
%it in a separate location. Remove any frames from the directory that you
%will not be analyzing downstream. 

%INPUT:
%dirname:Directory in which the image stack is saved.

%Calls upon:
%norm16bit.m

clear, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basename = '10012024_Exp2';
dirname=[basename '/' basename '_phase/' basename '_erased/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curdir=cd;

cd(dirname);
directory=dir('*.tif');
T=length(directory);
path(dirname,path)

%Pick the regions to erase  
%Change the third argument in the pickRegion function to the final
%pre-lysis frame
[rp1, rp2, count, im1]=pickRegion(dirname, directory, 11);
    
for t=1:T
    t
    
    for n=1:count
        %Load image
        imagename=directory(t).name;
        
        im=imread(imagename);
        [imM,imN]=size(im);
        
        [imcounts,bins]=hist(double(nonzeros(im1)));
        [~,mpos]=max(imcounts);
        val=bins(mpos);
        %maxval=bins(end-3);
        %val=maxval;
        im(rp1(n,2):rp2(n,2),rp1(n,1):rp2(n,1))=val*ones(rp2(n,2)-rp1(n,2)+1,rp2(n,1)-rp1(n,1)+1);
        delete(imagename);
        imwrite(im,imagename);
    end
end

close all

%% Functions
function [rp1, rp2, count, im1]=pickRegion(dirname, directory, t)
    
    %load image
    cd(dirname)
    imagename=directory(t).name;
    im1=imread(imagename);
    [imM,imN]=size(im1);

    ppix=0.5;
    im2=norm16bit(im1,ppix);

    %figure, imshowpair(imA, imB, 'montage')
           widthInches = 8;
       heightInches = 6;
       figure('Units', 'inches', 'Position', [1, 1, widthInches, heightInches])
    %figure,
    imshow(im1,stretchlim(im1)*65000)
    set(gcf,'Pointer','fullcross')
    hold on
    axis manual
    
    count=0;
    k=0;
    
    while k~=1
        count=count+1;
        k=waitforbuttonpress;
        point1=get(gca,'CurrentPoint');   
        finalRect=rbbox;                   
        %pause
        point2=get(gca,'CurrentPoint');    
        point1=point1(1,1:2);              
        point2=point2(1,1:2);
        point1(point1<1)=1;
        point2(point2<1)=1;
        
        if point1(2)>imM
            point1(2)=imM;
        end
        if point1(1)>imN
            point1(1)=imN;
        end
        if point2(2)>imM
            point2(2)=imM;
        end
        if point2(1)>imN
            point2(1)=imN;
        end
        p1=min(point1,point2);%Calculate locations
        p2=max(point1,point2);
        offset = abs(point1-point2);%And dimensions
        x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
        y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
        plot(x,y)

        rp1(count,:)=round(p1);
        rp2(count,:)=round(p2);
    end
 end
    
