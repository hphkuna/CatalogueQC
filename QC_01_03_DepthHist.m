%% --------
% QC_01_03_DepthHist(MeasCat)

% Displays histogram of event depths 

% Input: MeasCat - Measured catalogue


function QC_01_03_DepthHist(MeasCat, varargin)

PlotLimits = cell2mat(varargin);

close all

if ~isnumeric(MeasCat)
    error('Load catagoue in the correct format (see readme for description)')
end

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

figure('name', 'Depth Histogram', 'Position', [100, 100, 1049, 895])
DepthCats = 0:.5:ceil(max(EventDepth));
hist(EventDepth, DepthCats, 'Color', FirstColor)

if ~isempty(PlotLimits)
    xlim([PlotLimits(1) PlotLimits(2)])
else
    xlim([-.5 max(DepthCats)])
end

title('Depth Histogram', 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Depth [km]')
ylabel('Number of events')

h = findobj(gca,'Type','patch');
set(h, 'FaceColor', FirstColor)
set(h, 'EdgeColor', 'w');

print(gcf,'CurrentFigures/QC_01_03_DepthHistogram','-dpng', '-r300')



