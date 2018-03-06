%% --------
% QC_01_02_EpicPlotMagnitudeDepth(MeasCat)

% Displays simple epicentral plot with magnitude and depth of events

% Input: MeasCat - Measured catalogue


function QC_01_02_EpicPlotMagnitudeDepth(MeasCat, varargin)

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

% plot

if min(EventMag)<=0
    MagShift = abs(min(EventMag)) + .5;
else
    MagShift = 0;
end
DataScale = (EventMag + MagShift)*4;

figure('name', 'Epicentral Map w. magnitude and depth', 'Position', [100, 100, 1049, 895])
scatter(EventLon, EventLat, DataScale, EventDepth, 'filled', 'o')

if ~isempty(PlotLimits)
    ylim([PlotLimits(1) PlotLimits(3)])
    xlim([PlotLimits(2) PlotLimits(4)])
else
    xlim([min(EventLon) max(EventLon)])
    ylim([min(EventLat) max(EventLat)])
end

hold on
LonToLatRatio = (cosd((max(EventLat)-min(EventLat))/2)*111.3)/111.3;
daspect([1 LonToLatRatio 1])
h = colorbar;
ylabel(h, 'event depth')

title('Epicentral Map w. magnitude and depth', 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Longitude')
ylabel('Latitude')

print(gcf,'CurrentFigures/QC_01_02_EpicentralMapMagDepth','-dpng', '-r300')



