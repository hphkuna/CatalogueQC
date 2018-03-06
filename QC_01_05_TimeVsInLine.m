%% --------
% QC_01_05_TimeVsInLine(MeasCat)

% Displays Time vs InLine distance plot

% Input: MeasCat - Measured catalogue, coordinates of a rotation point and 
% the fault azimuth in degrees measured clockwise from north
% Rotation = [ZeroLat, ZeroLon, FaultAzimuth]


function QC_01_05_TimeVsInLine(MeasCat, Rotation, varargin)

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

%% rotation to in line direction

[EventDistance, EventAzimuth] = distance(Rotation(1), Rotation(2), EventLat, EventLon);

DegreeLength = 111.32;
EventDistance = EventDistance*DegreeLength;

XLineDistance = (sind(EventAzimuth-Rotation(3)).*EventDistance);
InLineDistance = sqrt(EventDistance.^2 - XLineDistance.^2);

%% display plot

% color definition
FirstColor = [0 .47 .95];
SecondColor = [.95 .47 0];
ThirdColor = [.33 .66 0];
Grey = [.7 .7 .7];

% plot

DataScale = zeros(length(EventMag),1) + 5;
DataScale(EventMag>=4) = 200;

% figure('name', 'Time vs. InLine Distance', 'Position', [100, 100, 1049, 895])
scatter(OriginTime, InLineDistance, DataScale, '.', 'MarkerEdgeColor', FirstColor)

if ~isempty(PlotLimits)
    xlim([min(OriginTime) max(OriginTime)])
    ylim([PlotLimits(1) PlotLimits(2)])
else
    xlim([min(OriginTime) max(OriginTime)])
    ylim([min(InLineDistance) max(InLineDistance)])
end

datetick('x',2)

xlabel('Date')
ylabel('In-line Distance [km]')
title('Time vs. InLine Distance', 'FontSize', 16, 'FontWeight', 'bold')


print(gcf,'CurrentFigures/QC_01_05_TimeVsInLine','-dpng', '-r300')




