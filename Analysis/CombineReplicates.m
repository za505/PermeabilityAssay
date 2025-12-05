%CombineReplicates
%Author: Zarina Akbary
%Combines data from nct.mat files

clear, close all

%INSTRUCTIONS FOR USE:
%run DataNormalization.m first

%INPUT
%basename: experiments of interest
%dirname: where .mat files are stored
%channels: list of directories containing fluorescent image stacks to quantify.

%OUTPUT:
%icell_intensity = cell x time matrix of fluoresc. intensities
%bg_intensity = 1 x time matrix of mean background intensities

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USER INPUT
basenames = {'10012024_Exp2', '10032024_Exp2', '10092024_Exp1'};
savename = ['tunicamycinSet_ponA_tunicamycin']; % name to save the .mat file
appendSet = 0; % no = 0, yes = 1

datadir = ['~/normalization']; %directory where nct.mat files are stored
save_plotdir = []; % directory to save plots
save_datadir = []; % directory to save .mat file

colors = lines(6);

% manually set the commonTime vector based on the time vector of the
% individual experiments
% commonTime = [0:59];
% cutoff = 60;

commonTime = [0:7, 12, 17];
cutoff = 10; % final time index 

if appendSet==1 %to add an experiment to an exisiting pooled dataset

    cd(save_datadir)
    load([savename '.mat'])
    n = length(data);

    cd(datadir)

    for i = 1:length(basenames)

        N = n + i;

        load([basenames{i} '_nct.mat'], 'fidx1', 'fidx2', 'time', 'icell_intensity', 'bg_intensity', 'pixels', 'B', 'lcell', 'new_lcell', 'adjintensity', 'correctedintensity', 'Cbl_exp', 'unb_frac', 'normintensity', 'normTime', 'slopes')
    
        % save the name of the experiment
        data(N).experiment = basenames{i};
    
        % filter for cells that do not have data at the 60 minute mark
        idx = find(sum(~isnan(normintensity), 2) > length(commonTime)/2);
        y = normintensity(idx, 1:cutoff);
    
        % save fluorescence values and calculate how many cells are in the
        % dataset
        combined_normTraces = [combined_normTraces; y];
    
        ntraces = length(idx);
        nsample = nsample + ntraces;
        data(N).ntraces = ntraces;
    
        % save the normTime and normintensity data
        data(N).idx = idx;
        data(N).normTime = normTime;
        data(N).normintensity = normintensity;
        %data(N).slopes = slopes;
        data(N).pixels = pixels;
        data(N).boundaries = B;
        data(N).new_lcell = new_lcell;
        data(N).final_prelysis_frame = fidx1;
        data(N).initial_postlysis_frame = fidx2;
        data(N).time = time;
        data(N).icell_intensity = icell_intensity;
        data(N).bg_intensity = bg_intensity;
        data(N).adjintensity = adjintensity;
        data(N).original_lcell = lcell;
        data(N).correctedintensity = correctedintensity;
        data(N).Cbl_exp = Cbl_exp;
        data(N).unbleached_fraction = unb_frac;
    
    end
    
else

    cd(datadir)
    data = struct;
    combined_normTraces = [];
    
    nsample = 0;

    for i = 1:length(basenames)

        load([basenames{i} '_nct.mat'], 'fidx1', 'fidx2', 'time', 'icell_intensity', 'bg_intensity', 'pixels', 'B', 'lcell', 'new_lcell', 'adjintensity', 'correctedintensity', 'Cbl_exp', 'unb_frac', 'normintensity', 'normTime', 'slopes')
    
        % save the name of the experiment
        data(i).experiment = basenames{i};
    
        % filter for cells that do not have data at the 60 minute mark
        idx = find(sum(~isnan(normintensity), 2) > length(commonTime)/2);
        y = normintensity(idx, 1:cutoff);
    
        % save fluorescence values and calculate how many cells are in the
        % dataset
        combined_normTraces = [combined_normTraces; y];
    
        ntraces = length(idx);
        nsample = nsample + ntraces;
        data(i).ntraces = ntraces;
    
        % save the normTime and normintensity data
        data(i).idx = idx;
        data(i).normTime = normTime;
        data(i).normintensity = normintensity;
        %data(i).slopes = slopes;
        data(i).pixels = pixels;
        data(i).boundaries = B;
        data(i).new_lcell = new_lcell;
        data(i).final_prelysis_frame = fidx1;
        data(i).initial_postlysis_frame = fidx2;
        data(i).time = time;
        data(i).icell_intensity = icell_intensity;
        data(i).bg_intensity = bg_intensity;
        data(i).adjintensity = adjintensity;
        data(i).original_lcell = lcell;
        data(i).correctedintensity = correctedintensity;
        data(i).Cbl_exp = Cbl_exp;
        data(i).unbleached_fraction = unb_frac;
    
    end

end

cd(save_datadir)
save([savename '.mat'], "data", "commonTime", "combined_normTraces", "nsample", "cutoff")

%% visualize all the traces from a given experiment
labels = {'rep1', 'rep2', 'rep3'};

for i = 1:length(data)
    subplot(1, 3 ,i)

    [nr, ~] = size(data(i).normintensity);

    for n = 1:nr
        plot(data(i).normTime, data(i).normintensity(n, :), 'Color', colors(i, :)), hold on
    end
    title(labels{i})
    ylabel('Normalized Fluorescence')
    xlabel('Time (minutes)')
    ylim([0 1.2])
end

