%FinalPlots.m
%Author: Zarina Akbary
%Purpose: to plot normalized fluorescence traces from pooled experimental
%replicates

%% mNG tunicamycin. Figure_3A. 
clear, close all

datadir = ['/Figure_3A'];
savedir = [];

basenames = {'tunicamycinSet_WT_untreated', 'tunicamycinSet_WT_tunicamycin', 'tunicamycinSet_ponA_untreated', 'tunicamycinSet_ponA_tunicamycin'}
labels = {'WT untreated', 'WT tunicamycin', 'ponA untreated', 'ponA tunicamycin'};

colors = lines(length(basenames));

% Define the width and height in inches
widthInches = 3 * 2.5;
heightInches = 2 * 2.5; % Adjust this value based on the desired aspect ratio

% first, plot the mean traces
figure('Units', 'inches', 'Position', [1, 1, widthInches, heightInches]), hold on
for i = 1:length(basenames)

    cd(datadir)
    load([basenames{i} '.mat'], 'commonTime', 'combined_normTraces', 'nsample')
            
    x = commonTime;
    y = combined_normTraces;
    ybar = mean(y, 1, 'omitnan');
    ysem = std(y, 0, 1, 'omitnan') ./ sqrt(nsample);
    y_upper = ybar + ysem;
    y_lower = ybar - ysem;

    plot(x, ybar, 'Color', colors(i, :), 'LineWidth', 1.2, 'DisplayName', labels{i})
    fill([x, fliplr(x)], [y_upper, fliplr(y_lower)], colors(i, :), 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', ' ')

end
legend('show', 'Location', 'best')
ylim([0 1.2])
xline(1, '--k', 'LineWidth', 1)
xline(5, '--k', 'LineWidth', 1)
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
xlim([0 60])

%cd(savedir)
%saveas(gcf, 'tunicamycinTreatment_means.svg')

% Define the width and height in inches
widthInches = 4.5 * 2.5;
heightInches = 1 * 2.5; % Adjust this value based on the desired aspect ratio

% Then, the individual data
for i = 1:length(basenames)

    cd(datadir)
    load([basenames{i} '.mat'], 'data')
       
    n1 = 1;
    n2 = length(data);

    figure('Units', 'inches', 'Position', [1, 1, widthInches, heightInches]), hold on
    
    for j = 1:length(data)

        x = data(j).normTime;
        y = data(j).normintensity;
        xtp = data(j).initial_postlysis_frame - data(j).final_prelysis_frame - 1;
        ybar = mean(y, 1, 'omitnan');
        ysem = std(y, 0, 1, 'omitnan') ./ data(j).ntraces;
        y_upper = ybar + ysem;
        y_lower = ybar - ysem;

        subplot(n1, n2, j), hold on
        plot(x, y, 'Color', colors(i, :), 'LineWidth', 1.2, 'DisplayName', ' ')

        plot(x, ybar, 'Color', 'black', 'LineWidth', 1.2, 'DisplayName', labels{i})
        fill([x, fliplr(x)], [y_upper, fliplr(y_lower)], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', ' ')
        
        title(['n = ' num2str(data(j).ntraces)])
        ylim([0 1.2])
        xline(1, '--k', 'LineWidth', 1)
        xline(xtp, '--k', 'LineWidth', 1)
        xlabel('Time (minutes)')
        ylabel('Normalized Fluorescence')
        xlim([0 60])
    end

    %cd(savedir)
    %saveas(gcf, [basenames{i} '_suplots.svg'])
 
end


%% GlpQ mNG set
clear, close all

datadir = ['/Figure_3B'];
savedir = [];

basenames = {'ER300_Buffer', 'ER300_GlpQ', 'ER608_Buffer', 'ER608_GlpQ'}
labels = {'WT Buffer', 'WT GlpQ', 'ponA Buffer', 'ponA GlpQ'};

colors = lines(length(basenames));

% Define the width and height in inches
widthInches = 3 * 2.5;
heightInches = 2 * 2.5; % Adjust this value based on the desired aspect ratio

% first, plot the mean traces
figure('Units', 'inches', 'Position', [1, 1, widthInches, heightInches]), hold on
for i = 1:length(basenames)

    cd(datadir)
    load([basenames{i} '.mat'], 'commonTime', 'combined_normTraces', 'nsample')
            
    x = commonTime;
    y = combined_normTraces;
    ybar = mean(y, 1, 'omitnan');
    ysem = std(y, 0, 1, 'omitnan') ./ sqrt(nsample);
    y_upper = ybar + ysem;
    y_lower = ybar - ysem;

    plot(x, ybar, 'Color', colors(i, :), 'LineWidth', 1.2, 'DisplayName', labels{i})
    fill([x, fliplr(x)], [y_upper, fliplr(y_lower)], colors(i, :), 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', ' ')

end
legend('show', 'Location', 'best')
ylim([0 1.2])
xline(1, '--k', 'LineWidth', 1)
xline(4.5, '--k', 'LineWidth', 1)
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
xlim([0 60])

%cd(savedir)
%saveas(gcf, 'GlpQ_mNG_means.svg')

% Then, the individual data
for i = 1:length(basenames)

    cd(datadir)
    load([basenames{i} '.mat'], 'data')
       
    n1 = 1;
    n2 = length(data);

    if i > 2
        % Define the width and height in inches
        widthInches = 4.5 * 2.5;
        heightInches = 1 * 2.5; % Adjust this value based on the desired aspect ratio
    else
        % Define the width and height in inches
        widthInches = 1.25 * 2.5;
        heightInches = 1 * 2.5; % Adjust this value based on the desired aspect ratio
    end


    figure('Units', 'inches', 'Position', [1, 1, widthInches, heightInches]), hold on
    
    for j = 1:length(data)

        x = data(j).normTime;
        y = data(j).normintensity;
        xtp = data(j).initial_postlysis_frame - data(j).final_prelysis_frame - 1;
        ybar = mean(y, 1, 'omitnan');
        ysem = std(y, 0, 1, 'omitnan') ./ data(j).ntraces;
        y_upper = ybar + ysem;
        y_lower = ybar - ysem;

        subplot(n1, n2, j), hold on
        plot(x, y, 'Color', colors(i, :), 'LineWidth', 1.2, 'DisplayName', ' ')

        plot(x, ybar, 'Color', 'black', 'LineWidth', 1.2, 'DisplayName', labels{i})
        fill([x, fliplr(x)], [y_upper, fliplr(y_lower)], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', ' ')
        
        title(['n = ' num2str(data(j).ntraces)])
        ylim([0 1.2])
        xline(1, '--k', 'LineWidth', 1)
        xline(xtp, '--k', 'LineWidth', 1)
        xlabel('Time (minutes)')
        ylabel('Normalized Fluorescence')
        xlim([0 60])
    end

    %cd(savedir)
    %saveas(gcf, [basenames{i} '_suplots.svg'])
 
end

%% GlpQ mUbq set
clear, close all

datadir = ['/Figure_3C'];
savedir = [];

basenames = {'ER607_Buffer', 'ER607_GlpQ'}
labels = {'ponA Buffer', 'ponA GlpQ'};

colors = lines(length(basenames));

% Define the width and height in inches
widthInches = 3 * 2.5;
heightInches = 2 * 2.5; % Adjust this value based on the desired aspect ratio

% first, plot the mean traces
figure('Units', 'inches', 'Position', [1, 1, widthInches, heightInches]), hold on
for i = 1:length(basenames)

    cd(datadir)
    load([basenames{i} '.mat'], 'commonTime', 'combined_normTraces', 'nsample')
            
    x = commonTime;
    y = combined_normTraces;
    ybar = mean(y, 1, 'omitnan');
    ysem = std(y, 0, 1, 'omitnan') ./ sqrt(nsample);
    y_upper = ybar + ysem;
    y_lower = ybar - ysem;

    plot(x, ybar, 'Color', colors(i, :), 'LineWidth', 1.2, 'DisplayName', labels{i})
    fill([x, fliplr(x)], [y_upper, fliplr(y_lower)], colors(i, :), 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', ' ')

end
legend('show', 'Location', 'best')
ylim([0 1.2])
xline(1, '--k', 'LineWidth', 1)
xline(4.5, '--k', 'LineWidth', 1)
xlabel('Time (minutes)')
ylabel('Normalized Fluorescence')
xlim([0 60])

%cd(savedir)
%saveas(gcf, 'GlpQ_mUbq_means.svg')


% Define the width and height in inches
widthInches = 4.5 * 2.5;
heightInches = 1 * 2.5; % Adjust this value based on the desired aspect ratio


% Then, the individual data
for i = 1:length(basenames)

    cd(datadir)
    load([basenames{i} '.mat'], 'data')
       
    n1 = 1;
    n2 = length(data);

    figure('Units', 'inches', 'Position', [1, 1, widthInches, heightInches]), hold on
    
    for j = 1:length(data)

        x = data(j).normTime;
        y = data(j).normintensity;
        xtp = data(j).initial_postlysis_frame - data(j).final_prelysis_frame - 1;
        ybar = mean(y, 1, 'omitnan');
        ysem = std(y, 0, 1, 'omitnan') ./ data(j).ntraces;
        y_upper = ybar + ysem;
        y_lower = ybar - ysem;

        subplot(n1, n2, j), hold on
        plot(x, y, 'Color', colors(i, :), 'LineWidth', 1.2, 'DisplayName', ' ')

        plot(x, ybar, 'Color', 'black', 'LineWidth', 1.2, 'DisplayName', labels{i})
        fill([x, fliplr(x)], [y_upper, fliplr(y_lower)], [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', ' ')
        
        title(['n = ' num2str(data(j).ntraces)])
        ylim([0 1.2])
        xline(1, '--k', 'LineWidth', 1)
        xline(xtp, '--k', 'LineWidth', 1)
        xlabel('Time (minutes)')
        ylabel('Normalized Fluorescence')
        xlim([0 60])
    end

    %cd(savedir)
    %saveas(gcf, [basenames{i} '_suplots.svg'])
 
end