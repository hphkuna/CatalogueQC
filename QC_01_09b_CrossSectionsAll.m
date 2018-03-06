%% --------
% QC_01_09_CrossSections(MeasCat)

% Displays sum of all sections with predefined parameters

% Input: MeasCat - Measured catalogue

function QC_01_09b_CrossSectionsAll(MeasCat, InLineDistance, SelectionDepth, Limits, varargin)

close all

if ~isnumeric(MeasCat)
    error('Load catagoue in the correct format (see readme for description)')
end

PlotLimits = cell2mat(varargin);

%% definition of variables

disp(' ')
warning('Parameters for the density plot might need some adjustment')
GridStepSection = 1; % distance between closest grid points in km
MaxDistSection = 1; % maximum distance to event in km
GridStepMap = .01; % distance between closest grid points in km
MaxDistMap = .1; % maximum distance to event in km

% definition of event origin time, latitude, longitude, depth and magnitude
OriginTime = datenum(MeasCat(:,1), MeasCat(:,2), MeasCat(:,3), MeasCat(:,4), MeasCat(:,5), MeasCat(:,6));
EventLat = MeasCat(:, 7);
EventLon = MeasCat(:, 8);
EventDepth = MeasCat(:, 9);
EventMag = MeasCat(:, 10);

% freeing up the memory
clear MeasCat

% degree length
LatDegreeLength = 111.132;
LonToLatRatio = (cosd(((max(EventLat)-min(EventLat))/2)+min(EventLat))*111.3)/111.3;
LonDegreeLength = LatDegreeLength*LonToLatRatio;

%% plot

% color definition
FirstColor = [0 .47 .95];
SecondColor = [.95 .47 0];
ThirdColor = [.33 .66 0];
Grey = [.7 .7 .7];

figure('name', 'Epicentral map with denoted cross sections', 'Position', [100, 100, 1049, 1500])
subplot(2,1,1)
plot(EventLon, EventLat, '.', 'Color', FirstColor, 'MarkerSize', 2)
if ~isempty(PlotLimits)
    xlim([PlotLimits(1) PlotLimits(3)])
    ylim([PlotLimits(5) PlotLimits(7)])
else
    xlim([min(EventLon) max(EventLon)])
    ylim([min(EventLat) max(EventLat)])
end
hold on
daspect([1 LonToLatRatio 1])
xlabel('Longitude')
ylabel('Latitude')

plot(Limits(1,:), Limits(2,:), '-', 'Color', SecondColor)
plot(Limits(3,:), Limits(4,:), '-', 'Color', SecondColor)
plot(Limits(5,:), Limits(6,:), '-', 'Color', SecondColor)
plot(Limits(7,:), Limits(8,:), '-', 'Color', SecondColor)
plot(Limits(9,:), Limits(10,:), '-', 'Color', SecondColor)

% display faults

subplot(2,1,2)
plot(InLineDistance, SelectionDepth, '.', 'MarkerSize', 2, 'Color', FirstColor, 'LineWidth', 2)

if ~isempty(PlotLimits)
    xlim([PlotLimits(2) PlotLimits(4)])
    ylim([PlotLimits(6) PlotLimits(8)])
else
    xlim([min(InLineDistance) max(InLineDistance)])
    ylim([0 max(EventDepth)])
end

set(gca, 'YDir', 'reverse')
daspect([1 1 1])

xlabel('In-line Distance [km]')
ylabel('Hypocenter Depth [km]')
% title('In-line distance vs. depth', 'FontSize', 16, 'FontWeight', 'bold')

print(gcf,'CurrentFigures/QC_01_09_InLineVsDepth_section','-dpng', '-r300')

%% event density

%% epicental map density

[Density, X, Y] = xDensityPlot(EventLon, EventLat, GridStepMap, MaxDistMap);

% plot
figure('name', 'Event Density Map', 'Position', [100, 100, 1049, 1500])

subplot(2,1,1)
surf(X, Y, Density, 'LineStyle', 'none', 'FaceColor', 'interp')
daspect([1 LonToLatRatio 1])
view(0, 90)
h = colorbar;
% set(h, 'Position', [.92 .6 .02 .2])
ylabel(h, 'log10 of Number of Events')
grid off

if ~isempty(PlotLimits)
    xlim([PlotLimits(1) PlotLimits(3)])
    ylim([PlotLimits(5) PlotLimits(7)])
else
    xlim([min(EventLon) max(EventLon)])
    ylim([min(EventLat) max(EventLat)])
end

hold on
plot3(Limits(1,:), Limits(2,:), ones(5, 1) + 5, '-', 'Color', SecondColor)
plot3(Limits(3,:), Limits(4,:), ones(5, 1) + 5, '-', 'Color', SecondColor)
plot3(Limits(5,:), Limits(6,:), ones(5, 1) + 5, '-', 'Color', SecondColor)
plot3(Limits(7,:), Limits(8,:), ones(5, 1) + 5, '-', 'Color', SecondColor)
plot3(Limits(9,:), Limits(10,:), ones(5, 1) + 5, '-', 'Color', SecondColor)

disp(['Circle radius: ' num2str(MaxDistSection)])
disp(['Grid step: ' num2str(GridStepSection)])
disp(' ')

% title('Event Density Map', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Longitude')
xlabel('Latitude')
    

%% cross section

[Density, X, Y] = xDensityPlot(InLineDistance, SelectionDepth, GridStepSection, MaxDistSection);   

% plot
subplot(2,1,2)
surf(X, Y, Density, 'LineStyle', 'none', 'FaceColor', 'interp')
view(0,-90)
daspect([1 1 1])
h = colorbar;
% set(h, 'Position', [.92 .1 .02 .2])
ylabel(h, 'log10 of Number of Events')

if ~isempty(PlotLimits)
    xlim([PlotLimits(2) PlotLimits(4)])
    ylim([PlotLimits(6) PlotLimits(8)])
else
    xlim([min(InLineDistance) max(InLineDistance)])
    ylim([0 max(EventDepth)])
end

disp(['Circle radius: ' num2str(MaxDistSection)])
disp(['Grid step: ' num2str(GridStepSection)])
disp(' ')
grid off

% title('Event Density Map', 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Hypocenter Depth [km]')
ylabel('In-line Distance [km]')

print(gcf,'CurrentFigures/QC_01_09b_InLineVsDepthDens_sectionAll','-dpng', '-r300')



















