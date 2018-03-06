%% --------
% QC_01_06_TimeHistogram(MeasCat)

% Time Histogram

% Input: MeasCat - Measured catalogue


function QC_01_06_TimeHistogram(MeasCat, varargin)

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

%% plot

% color definition
FirstColor = [0 .47 .95];
SecondColor = [.95 .47 0];
ThirdColor = [.33 .66 0];
Grey = [.7 .7 .7];

% plot

figure('name', 'Time Histogram', 'Position', [100, 100, 1049, 895])
DateBins = floor(min(OriginTime)):.1:ceil(max(OriginTime));
hist(OriginTime, DateBins, 'Color', FirstColor)
datetick('x', 2)
xlim([floor(min(OriginTime))-5 ceil(max(OriginTime))])
title('Time Histogram', 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Date')
ylabel('Number of events')

h = findobj(gca,'Type','patch');
set(h, 'FaceColor', FirstColor)
% set(h, 'EdgeColor', 'w');

if ~isempty(PlotLimits)
    xlim([min(OriginTime) max(OriginTime)])
    ylim([PlotLimits(1) PlotLimits(2)])
end

print(gcf,'CurrentFigures/QC_01_06_TimeHistogram','-dpng', '-r300')









