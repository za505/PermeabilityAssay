%ImageAlign.m
%Rico Rojas
%Registers image stacks using a common reference stack.  
%Calls on dftregistration.m
%
%INSTRUCTIONS FOR USE:
%For a given movie, split the color channels and save each image stack you 
%want to register (the phase and fluorescence channel) in separate 
%directories. Also, crop a portion of the phase movie (which contains 
%features from which it is possible to track drift) and save the reference 
%image stack in its own directory. The program will create new directories 
%for each of the image stacks to be registered in the same
%parent directory. The reference stacks and the stacks to be registered
%should have the same number of images.
%
%INPUT
%basename: choose a name for the directories in which the aligned images
%          will be saved.
%dirname:  cell array containing the paths of the directories in which the
%          image stacks to be registered are.
%regname:  path of reference directory.
%
%CALLS ON:
%dftregistration.m
%imtranslate.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%User Input
basename='10012024_Exp2'
dirname={[basename '/' basename '_phase/' basename '_full']; [basename '/' basename '_mNeonGreen/' basename '_full']};
regname=[basename '/' basename '_phase/' basename '_crop']; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

workdir=pwd;
ld=length(dirname); %number of channels being aligned

cd(regname)
directory=dir('*.tif');
path(regname,path)
T=length(directory); %number of frames in each image stack
tf=T*ld; %number of total frames to be aligned

imagename=directory(1).name; %load the first frame in the register image stack
RefIm=imread([regname '/' imagename]);

count=0; %to keep track of how many frames have been aligned
for i=1:ld
    
    cd(dirname{i})
    path(dirname{i},path)
    directory2=dir('*.tif');
    cd('../')
    
    mkdir([basename '_aligned'])
    cd(['./' basename '_aligned'])
    
    shft=zeros(T,2);
    
    for t=1:T
        count=count+1;
        pd=count/tf*100;
        pds=sprintf(['%4.1f'],pd);
        [pds '%'] 

        imagename=directory(t).name;
        TestIm=imread([regname '/' imagename]);

        [output NewImFT]=dftregistration(fft2(RefIm),fft2(TestIm), 1); %all subsequent images are aligned to the first one in the stack
        NewIm=abs(ifft2(NewImFT));

        shft(t+1,1)=output(3);
        shft(t+1,2)=output(4);

        imagename=directory2(t).name;

        I=imread([dirname{i} '/' imagename]);
        [counts,bins]=imhist(I);
        [~,maxpos]=max(counts);
        padcolor=bins(maxpos);

        I=imtranslate(I,shft(t+1,:),padcolor);

        b=sprintf(['%4.4d'],t);  
        savename=[basename '_a' b '.tif'];

        imwrite(I,savename);
    end
end



