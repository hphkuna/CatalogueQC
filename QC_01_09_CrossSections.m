%% --------
% QC_01_09_CrossSections(MeasCat)

% Displays cross sections

% Input: MeasCat - Measured catalogue, coordinates of a rotation point and 
% the fault azimuth in degrees measured clockwise from north


function [InLineDistance, EventDepth, Limits] = QC_01_09_CrossSections(MeasCat, number, par, PlotOut, varargin)

% close all

if ~isnumeric(MeasCat)
    error('Load catagoue in the correct format (see readme for description)')
end

PlotLimits = cell2mat(varargin);

%% definition of variables

disp(' ')
warning('Parameters for the density plot might need some adjustment')
GridStepSection = 1; % distance between closest grid points in km
MaxDistSection = 1; % maximum distance to event in km
GridStepMap = .025; % distance between closest grid points in km
MaxDistMap = .1; % maximum distance to event in km

disp(' ')
disp('Magnitude range can be constrained')
disp(' ')
ConstrainMagnitude = -5;

LowerLeftCor = [par(1) par(2)]; % latitude and longitude
RotationAngle = par(3); % in degrees
SectionWidth = par(4);
SectionLength = 150; % in km

% definition of event origin time, latitude, longitude, depth and magnitude
OriginTime = datenum(MeasCat(:,1), MeasCat(:,2), MeasCat(:,3), MeasCat(:,4), MeasCat(:,5), MeasCat(:,6));
EventLat = MeasCat(:, 7);
EventLon = MeasCat(:, 8);
EventDepth = MeasCat(:, 9);
EventMag = MeasCat(:, 10);
StrongEv = find(EventMag>=4);

% freeing up the memory
clear MeasCat

% constrain with magnitude??
EventLat = EventLat(EventMag>=ConstrainMagnitude);
EventLon = EventLon(EventMag>=ConstrainMagnitude);
EventDepth = EventDepth(EventMag>=ConstrainMagnitude);

% degree length
LatDegreeLength = 111.132;
LonToLatRatio = (cosd(((max(EventLat)-min(EventLat))/2)+min(EventLat))*111.3)/111.3;
LonDegreeLength = LatDegreeLength*LonToLatRatio;

%% rotation to in line direction

EventLatKm = EventLat*LatDegreeLength;
EventLonKm = EventLon*LonDegreeLength;

LowerLeftCorLatKm = LowerLeftCor(1)*LatDegreeLength;
LowerLeftCorLonKm = LowerLeftCor(2)*LonDegreeLength;

Diff1 = [EventLatKm - LowerLeftCorLatKm, EventLonKm - LowerLeftCorLonKm]';

RotationMatrix = [cosd(-RotationAngle) -sind(-RotationAngle); sind(-RotationAngle) cosd(-RotationAngle)];

Diff2 = RotationMatrix * Diff1;

InLineDistance = Diff2(1, :);
XLineDistance = Diff2(2, :);

% select events in cross section
IndexInSection = find(XLineDistance<=SectionWidth & XLineDistance>=0);

%% create mesh

UpperLeftCorLatKm = LowerLeftCorLatKm + SectionLength;
UpperLeftCorLonKm = LowerLeftCorLonKm;

LowerRightCorLatKm = LowerLeftCorLatKm;
LowerRightCorLonKm = LowerLeftCorLonKm + SectionWidth;

UpperRightCorLatKm = UpperLeftCorLatKm;
UpperRightCorLonKm = LowerRightCorLonKm;

BoundLat = [LowerLeftCorLatKm LowerRightCorLatKm UpperRightCorLatKm UpperLeftCorLatKm LowerLeftCorLatKm];
BoundLon = [LowerLeftCorLonKm LowerRightCorLonKm UpperRightCorLonKm UpperLeftCorLonKm LowerLeftCorLonKm];

% rotate mesh
Diff1 = [BoundLat - LowerLeftCorLatKm; BoundLon - LowerLeftCorLonKm];

RotationMatrix = [cosd(RotationAngle) -sind(RotationAngle); sind(RotationAngle) cosd(RotationAngle)];

Diff2 = RotationMatrix * Diff1;

BoundLatKm = BoundLat(1) + Diff2(1, :);
BoundLonKm = BoundLon(1) + Diff2(2, :);

BoundLat = BoundLatKm / LatDegreeLength;
BoundLon = BoundLonKm / LonDegreeLength;
Limits = [BoundLon; BoundLat];

InLineDistance = InLineDistance*cosd(-RotationAngle); % label with km towards North
InLineDistance = InLineDistance(IndexInSection)';
EventDepthSel = EventDepth(IndexInSection);

%% plot

if PlotOut == 1

% color definition
FirstColor = [0 .47 .95];
SecondColor = [.95 .47 0];
ThirdColor = [.33 .66 0];
Grey = [.7 .7 .7];

figure('name', 'Epicentral map with denoted cross sections', 'Position', [100, 100, 1049, 1500])
subplot(2,1,1)


plot(EventLon(EventDepth>10), EventLat(EventDepth>10), 'ro', 'MarkerSize', 1)
hold on
plot(EventLon(EventDepth<10), EventLat(EventDepth<10), '.', 'Color', FirstColor, 'MarkerSize', 1)
% xlim([min([EventLon; BoundLon']) max([EventLon; BoundLon'])])
% ylim([min([EventLat; BoundLat']) max([EventLat; BoundLat'])])

if ~isempty(PlotLimits)
    xlim([PlotLimits(1) PlotLimits(3)])
    ylim([PlotLimits(5) PlotLimits(7)])
else
    xlim([min(EventLon) max(EventLon)])
    ylim([min(EventLat) max(EventLat)])
end

hold on
daspect([1 LonToLatRatio 1])
% title('Basic Epicentral Map', 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Longitude')
ylabel('Latitude')

plot(BoundLon, BoundLat, '-', 'Color', SecondColor)
plot(EventLon(IndexInSection), EventLat(IndexInSection), 'o', 'MarkerSize', .1, 'Color', SecondColor)

subplot(2,1,2)
plot(InLineDistance, EventDepthSel, '.', 'MarkerSize', 5, 'Color', FirstColor)

if ~isempty(PlotLimits)
    xlim([PlotLimits(2) PlotLimits(4)])
    ylim([PlotLimits(6) PlotLimits(8)])
else
    xlim([min(InLineDistance) max(InLineDistance)])
    ylim([0 max(EventDepthSel)])
end

set(gca, 'YDir', 'reverse')
daspect([1 1 1])

xlabel('In-line Distance [km]')
ylabel('Hypocenter Depth [km]')
% title('In-line distance vs. depth', 'FontSize', 16, 'FontWeight', 'bold')

print(gcf,['CurrentFigures/QC_01_09_InLineVsDepth_section' num2str(number)],'-depsc', '-r300')


%% density of events

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

BoundZ = ones(length(BoundLon), 1) + 5;

hold on
plot3(BoundLon, BoundLat, BoundZ, '-', 'Color', SecondColor)

disp(['Circle radius: ' num2str(MaxDistMap)])
disp(['Grid step: ' num2str(GridStepMap)])
disp(' ')

% title('Event Density Map', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Longitude')
xlabel('Latitude')
    

%% cross section

[Density, X, Y] = xDensityPlot(InLineDistance, EventDepthSel, GridStepSection, MaxDistSection);   

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
    ylim([0 max(EventDepthSel)])
end

disp(['Circle radius: ' num2str(MaxDistSection)])
disp(['Grid step: ' num2str(GridStepSection)])
disp(' ')
grid off

% title('Event Density Map', 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Hypocenter Depth [km]')
ylabel('In-line Distance [km]')

print(gcf,['CurrentFigures/QC_01_09_InLineVsDepth_sectionD' num2str(number)],'-dpng', '-r300')
    

end





