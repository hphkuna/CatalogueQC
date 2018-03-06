%% --------
% QC_01_10_TimeVsDepth(MeasCat)

% Displays scatter plot of tive vs depth

% Input: MeasCat - Measured catalogue


function QC_01_10_TimeVsDepth(MeasCat, varargin)

close all

if ~isnumeric(MeasCat)
    error('Load catagoue in the correct format (see readme for description)')
end

PlotLimits = cell2mat(varargin);

%% definition of variables

% definition of event origin time, latitude, longitude, depth and magnitude
OriginTime = datenum(MeasCat(:,1), MeasCat(:,2), MeasCat(:,3), MeasCat(:,4), MeasCat(:,5), MeasCat(:,6));
EventLat = MeasCat(:, 7);
EventLon = MeasCat(:, 8);
EventDepth = MeasCat(:, 9);
EventMag = MeasCat(:, 10);

% freeing up the memory
clear MeasCat

%% display plot

% color definition
FirstColor = [0 .47 .95];
SecondColor = [.95 .47 0];
ThirdColor = [.33 .66 0];
Grey = [.7 .7 .7];

% plot
figure('name', 'Time Vs. Depth diagram', 'Position', [100, 100, 1049, 895])
plot(OriginTime, EventDepth, '.', 'Color', FirstColor, 'MarkerSize', 4)

if ~isempty(PlotLimits)
    ylim([PlotLimits(1) PlotLimits(3)])
    xlim([PlotLimits(2) PlotLimits(4)])
else
    xlim([min(OriginTime) max(OriginTime)])
    ylim([0 max(EventDepth)])
end

title('Time Vs. Depth diagram', 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Time')
ylabel('Depth')
datetick('x',2)
set(gca, 'YDir', 'reverse')

print(gcf,'CurrentFigures/QC_01_10_TimeVsDepth','-dpng', '-r300')











