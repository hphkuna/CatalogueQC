%% --------
% QC_01_01_EpicentralPlot(MeasCat)

% Displays simple epicentral plot of events in a catalogue and density plot
% of events

% Input: MeasCat - Measured catalogue


function QC_01_01_EpicentralPlot(MeasCat, varargin)

close all

if ~isnumeric(MeasCat)
    error('Load catagoue in the correct format (see readme for description)')
end

PlotLimits = cell2mat(varargin);

%% definition of variables for density plot

disp(' ')
warning('Parameters for the density plot might need some adjustment')
GridStep = .1; % distance between closest grid points in degrees
MaxDist = .2; % maximum distance to event in degrees

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
figure('name', 'Basic Epicentral Map', 'Position', [100, 100, 1049, 895])
plot(EventLon, EventLat, '.', 'Color', FirstColor, 'MarkerSize', 4)

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

title('Basic Epicentral Map', 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Longitude')
ylabel('Latitude')

print(gcf,'CurrentFigures/QC_01_01_EpicentralMap','-dpng', '-r300')


%% density of events

[Density, X, Y] = xDensityPlot(EventLon, EventLat, GridStep, MaxDist);

% plot
figure('name', 'Event Density Map', 'Position', [100, 100, 1049, 1500])

surf(X, Y, Density, 'LineStyle', 'none', 'FaceColor', 'interp')
daspect([1 LonToLatRatio 1])
view(0, 90)
h = colorbar;
% set(h, 'Position', [.92 .6 .02 .2])
ylabel(h, 'log10 of Number of Events')
grid off

if ~isempty(PlotLimits)
    ylim([PlotLimits(1) PlotLimits(3)])
    xlim([PlotLimits(2) PlotLimits(4)])
else
    xlim([min(EventLon) max(EventLon)])
    ylim([min(EventLat) max(EventLat)])
end

disp(['Circle radius: ' num2str(MaxDist)])
disp(['Grid step: ' num2str(GridStep)])
disp(' ')

title('Event Density Map', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Longitude')
xlabel('Latitude')

print(gcf,'CurrentFigures/QC_01_01_EpicentralMapDensity','-dpng', '-r300')


%% counts and basic catalogue info

NumberOfEvents = length(OriginTime);
DateMin = datestr(min(OriginTime));
DateMax = datestr(max(OriginTime));
LatMin = min(EventLat);
LatMax = max(EventLat);
LonMin = min(EventLon);
LonMax = max(EventLon);
MinDepth = min(EventDepth);
MaxDepth = max(EventDepth);

%% print

disp('----------')
disp(' ')
disp('QC_01_01_EpicentralPlot - Output:')
disp(' ')
disp(['Number of events: ' num2str(NumberOfEvents)])
disp(['Start date: ' num2str(DateMin)])
disp(['End date: ' num2str(DateMax)])
disp(['Minimum latitude: ' num2str(LatMin)])
disp(['Maximum latitude: ' num2str(LatMax)])
disp(['Minimum longitude: ' num2str(LonMin)])
disp(['Maximum longitude: ' num2str(LonMax)])
disp(['Minimum depth: ' num2str(MinDepth)])
disp(['Maximum depth: ' num2str(MaxDepth)])
disp(' ')
disp('----------')







