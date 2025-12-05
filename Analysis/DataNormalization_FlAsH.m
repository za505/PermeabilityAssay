%DataNormalization_FlAsH.m
%Author: Zarina Akbary
%Corrects for photobleaching and normalizes data to prel-lysis frame

clear, close all

%INSTRUCTIONS FOR USE:
%run DiffusionMeasure.m first

%INPUT
%basename: experiments of interest
%datadir: where .mat input files are stored
%savedir: where .mat output files are stored
%autofluorfile: .mat with uninduced fluorescence values

%OUTPUT:
%normintensity = photobleach corrected, background subtracted fluor cell x
%time fluor matrix normalized to final prelysis value
%normTime = 1 x time matrix from final prelysis frame onwards

%% Save uninduced autofluorescence as it's own .mat file
% datadir = ['~/fluorescence'];
% cd(datadir)
% basename = '08132025_Exp3';
% 
% load([basename '_dm.mat'], 'icell_intensity', 'bg_intensity', 'time')
% adjintensity = icell_intensity - bg_intensity;
% mean_uninduced = mean(adjintensity, 1, 'omitnan');
% autoTime = time;
% 
% save('uninduced_FlAsH_noWash.mat')

%% User Inputs
datadir = ['~/fluorescence'];
savedir = ['~/normalization'];
%autofluorfile = ['uninduced_FlAsH_slrR.mat'];
%autofluorfile = ['uninduced_FlAsH_noWash.mat'];
%autofluorfile = ['uninduced_FlAsH_ponA_Exp2.mat'];
autofluorfile = ['uninduced_FlAsH.mat'];

% set rho based on fluorophore
rho_FlAsH = 0.02;

fidx1 = 5; %final pre-lysis frame
fidx2 = 8; %initial post-lysis frame

basename = '10012024_Exp2'; %Name used to save file.

% time vector for autofluor (if it doesn't have one saved in .mat file) 
autoTime = [[0:11], [16:5:72]];
%autoTime = [[0:11], [16:5:21]];

%% Perform Normalization and Correction
cd(datadir)

load([basename '_dm.mat'], 'icell_intensity', 'time', 'bg_intensity', 'lcell', 'acell', 'wcell', 'pixels', 'B')

 % interpolate post-lysis length values
[nr, ~] = size(lcell);
new_lcell = lcell;
    
for n = 1:nr
    
    y = lcell(n, fidx2:end);
    x = find(~isnan(y));
    xq = find(isnan(y));

    v = y(x);
    vq = interp1(x, v, xq);

    y(xq) = vq;
    new_lcell(n, fidx2:end) = y;
end

% interpolate autofluorescence value (nonspecific binding)
load(autofluorfile, 'mean_uninduced');
[~, nc] = size(mean_uninduced);
autofluorescence = interp1(autoTime, mean_uninduced, time);
%autofluorescence = mean_uninduced;

% subtract the background intensity and autofluorescence
adjintensity = icell_intensity - bg_intensity;
adjintensity = adjintensity - autofluorescence;

% set negative values to zero (the background intensity in the reference
% region was higher than the background of the area where the cell was
adjintensity(adjintensity < 0) = 0;

% normalize by pre-lysis fluorescence
norm_intensity = adjintensity ./ adjintensity(:, fidx1);

% interpolate the lysis frames
[nrow, ncol] = size(norm_intensity);
x = setdiff(1:ncol, [fidx1+1:fidx2-1]);
xq = time(fidx1+1:fidx2-1);

for n = 1:nrow
    v = norm_intensity(n, x);
    vq = interp1(x,v,xq);
    norm_intensity(n, fidx1+1:fidx2-1) = vq;
end

% perform photobleach correction on post-lysis frames
y = norm_intensity(:, fidx2:end);
[nr, nc] = size(y);

% preallocate array
Cnew = nan(nr, nc);
Cnew(:, 1) = y(:, 1);

Cbl_exp = nan(nr, nc);
Cbl_exp(:, 1) = 0; 

unb_frac = nan(nr, nc);
unb_frac(:, 1) = 1; 

for n = 1:nr
    for f = 1:nc-1

        dCB = y(n,f) * rho_FlAsH;
        dCT = y(n, f+1) - y(n, f);

        if dCB > abs(dCT) %to prevent over-correction (sometimes dCB is only a bit bigger than abs(dCT), which is basically abs(dCT)==dCB anyways. In cases where the discrepancy is bigger, I don't know what's going on, but I don't think it's permeability) 
           dCP = 0;
           Cbl_exp(n,f+1)=Cbl_exp(n,f)+abs(dCT);
        else
           dCP = dCT + dCB;
           Cbl_exp(n,f+1)=Cbl_exp(n,f)+dCB;
        end
        
        if dCP > 0 
            dCP = 0;
        end

        Cnew(n, f+1) = Cnew(n, f) + dCP;
        
        unb_frac(n,f+1)=(y(n,f+1))/(y(n,f+1)+Cbl_exp(n,f+1));%Calculate the new fraction of unbleached fluorophores

    end
end

correctedintensity = norm_intensity;
correctedintensity(:, fidx2:end) = Cnew;

% subset the normalized intensity to include only the final prelysis frame
% onwards
normintensity = correctedintensity(:, fidx1:end);
normTime = time(fidx1:end) - time(fidx1);

figure(1)
plot(normTime, normintensity, '-r')
ylim([0 1.2])
xline(1, '--k')
xline(3, '--k')
xlabel('Time (minutes)')
ylabel('Normalzied Fluorescence (AU)')

cd(savedir)
save([basename '_nct.mat'])

