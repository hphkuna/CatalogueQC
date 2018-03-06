%% --------
% QC_01_04_InLineVsDepth(MeasCat)

% Displays cross section of the fault

% Input: MeasCat - Measured catalogue, coordinates of a rotation point and 
% the fault azimuth in degrees measured clockwise from north
% Rotation = [ZeroLat, ZeroLon, FaultAzimuth]


function QC_01_04_InLineVsDepth(MeasCat, Rotation, varargin)

PlotLimits = cell2mat(varargin);

% close all

if ~isnumeric(MeasCat)
    error('Load catagoue in the correct format (see readme for description)')
end


%% definition of variables

disp(' ')
warning('Parameters for the density plot might need some adjustment')
GridStep = .5; % distance between closest grid points in km
MaxDist = 3; % maximum distance to event in km

disp(' ')
disp('Magnitude range can be constrained')
disp(' ')
ConstrainMagnitude = -5;

%% definition of variables

% definition of event origin time, latitude, longitude, depth and magnitude
OriginTime = datenum(MeasCat(:,1), MeasCat(:,2), MeasCat(:,3), MeasCat(:,4), MeasCat(:,5), MeasCat(:,6));
EventLat = MeasCat(:, 7);
EventLon = MeasCat(:, 8);
EventDepth = MeasCat(:, 9);
EventMag = MeasCat(:, 10);

% freeing up the memory
clear MeasCat

% constrain with magnitude??
EventLat = EventLat(EventMag>=ConstrainMagnitude);
EventLon = EventLon(EventMag>=ConstrainMagnitude);
EventDepth = EventDepth(EventMag>=ConstrainMagnitude);

%% rotation to in line direction

[EventDistance, EventAzimuth] = distance(Rotation(1), Rotation(2), EventLat, EventLon);

DegreeLength = 111.32;
EventDistance = EventDistance*DegreeLength;

XLineDistance = (sind(EventAzimuth-Rotation(3)).*EventDistance);
InLineDistance = sqrt(EventDistance.^2 - XLineDistance.^2);

%% plot

% color definition
FirstColor = [0 .47 .95];
SecondColor = [.95 .47 0];
ThirdColor = [.33 .66 0];
Grey = [.7 .7 .7];

figure('name', 'In-line distance vs. depth', 'Position', [100, 100, 1049, 895])
plot(InLineDistance, EventDepth, '.', 'MarkerSize', 8, 'Color', FirstColor)
set(gca, 'YDir', 'reverse')
daspect([3 1 1])

a = [InLineDistance, EventDepth];
save('Depth_Inline', 'a', '-ascii')
a = [XLineDistance, EventDepth];
save('Depth_Xline', 'a', '-ascii')

if ~isempty(PlotLimits)
    xlim([PlotLimits(1) PlotLimits(3)])
    ylim([PlotLimits(2) PlotLimits(4)])
else
    xlim([min(InLineDistance) max(InLineDistance)])
    ylim([0 max(EventDepth)])
end

xlabel('In-line Distance [km]')
ylabel('Hypocenter Depth [km]')
title('In-line distance vs. depth', 'FontSize', 16, 'FontWeight', 'bold')

print(gcf,'CurrentFigures/QC_01_04_InLineVsDepth','-dpng', '-r300')

%% density of events

[Density, X, Y] = xDensityPlot(InLineDistance, EventDepth, GridStep, MaxDist); 

% plot
figure('name', 'Event Density Map', 'Position', [100, 100, 1049, 895])

surf(X, Y, Density, 'LineStyle', 'none', 'FaceColor', 'interp')
view(0, -90)
daspect([1 1 1])
h = colorbar;
set(h, 'Position', [.92 .3 .02 .4])
ylabel(h, 'log10 of Number of Events')
daspect([3 1 1])

if ~isempty(PlotLimits)
    xlim([PlotLimits(1) PlotLimits(3)])
    ylim([PlotLimits(2) PlotLimits(4)])
else
    xlim([min(InLineDistance) max(InLineDistance)])
    ylim([0 max(EventDepth)])
end

disp(['Circle radius: ' num2str(MaxDist)])
disp(['Grid step: ' num2str(GridStep)])
disp(' ')
grid off

title('Event Density Map', 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Hypocenter Depth [km]')
ylabel('In-line Distance [km]')

print(gcf,'CurrentFigures/QC_01_04_InLineVsDepthDensity','-dpng', '-r300')

