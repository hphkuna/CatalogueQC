%% --------
% QC_03_03_McSpaceVariation(MeasCat)

% Computes Mc spatial variation by MAXC technique

% Input: MeasCat - Measured catalogue


function QC_03_03_McSpaceVariation(MeasCat, varargin)

close all

if ~isnumeric(MeasCat)
    error('Load catagoue in the correct format (see readme for description)')
end

PlotLimits = cell2mat(varargin);

%% definition of variables for density plot

disp(' ')
warning('Parameters for the density plot might need some adjustment')
GridStep = .02; % distance between closest grid points in degrees
MaxDist = .1; % maximum distance to event in degrees
MinNumberOfEvents = 40; % minimum number of events for Magnitude completeness computation

%% definition of variables

% definition of event origin time, latitude, longitude, depth and magnitude
OriginTime = datenum(MeasCat(:,1), MeasCat(:,2), MeasCat(:,3), MeasCat(:,4), MeasCat(:,5), MeasCat(:,6));
EventLat = MeasCat(:, 7);
EventLon = MeasCat(:, 8);
EventDepth = MeasCat(:, 9);
EventMag = MeasCat(:, 10);

% freeing up the memory
clear MeasCat

%% time variation of Mc

% definition of MagBins
MinMag = (floor(min(EventMag)*10))/10;
MaxMag = (ceil(max(EventMag)*10))/10;
MagBins = MinMag:.1:MaxMag;

% definition of vectors for grid mesh
LonToLatRatio = (cosd((max(EventLat)-min(EventLat))/2)*111.3)/111.3;

LatGrid = (floor(min(EventLat)):GridStep:ceil(max(EventLat)))';
LonGrid = (floor(min(EventLon)):GridStep*1/LonToLatRatio:ceil(max(EventLon)))';

LatVector = sort(repmat(LatGrid, length(LonGrid),1));
LonVector = repmat(LonGrid, length(LatGrid),1);

% distance function
EvalFct = @(LatGrid, LonGrid, EventLat, EventLon) ...
    sqrt((EventLat-LatGrid).^2 + (EventLon-LonGrid).^2);

% grid mesh
[LatMesh1, LatMesh2] = meshgrid(LatVector, EventLat);
[LonMesh1, LonMesh2] = meshgrid(LonVector, EventLon);

% evaluation of the distance function
Distances = EvalFct(LatMesh1, LonMesh1, LatMesh2, LonMesh2);

% replicate vector of magnitudes
EventMagMatrix = repmat(EventMag, 1, length(LatVector));

% select and compute histogram of magnitudes belonging to events closer
% that distance MaxDist
Indexes = (Distances<MaxDist');
Histograms = hist(EventMagMatrix.*Indexes, MagBins, 1);

% There is many zero entries that get summed into histogram and the
% histogram must be corrected for them. If the catalugue contains zero
% magnitudes, the zeros are subtracted from the zero magnitude bin,
% otherwise from the lowest magnitude bin
Correction2Zeros = length(EventMag) - sum(Indexes,1);

% Just making sure that there are no error accumulated over steps
MagBins = round(MagBins*10)/10;
if MinMag<=0
    Histograms(MagBins==0,:) = Histograms(MagBins==0,:) - Correction2Zeros;
else
    Histograms(1,:) = Histograms(1,:) - Correction2Zeros;
end

% Indexes of maximums in histograms
[~, MaxIndex] = max(Histograms);
% And respective magnitude bins
MagnitudeCompleteness = MagBins(MaxIndex);

% Number of events at each gridpoint. Only those with number greater than
% MinNumberOfEvents are considered
NumberOfEvents = sum(Histograms);
NumberOfEvents(NumberOfEvents<MinNumberOfEvents) = 0;
NumberOfEvents(NumberOfEvents>=MinNumberOfEvents) = 1;

% Consider only grid points with number of events higher than MinNumberOfEvents
MagnitudeCompleteness = MagnitudeCompleteness.*NumberOfEvents;

% Reshape the final matrix
MagnitudeCompleteness = reshape(MagnitudeCompleteness, length(LonGrid), length(LatGrid));
MagnitudeCompleteness(MagnitudeCompleteness==0) = NaN;


%% display plot

% color definition
FirstColor = [0 .47 .95];
SecondColor = [.95 .47 0];
ThirdColor = [.33 .66 0];
Grey = [.7 .7 .7];

% plot
figure('name', 'Event Density Map', 'Position', [100, 100, 1049, 895])

surf(LatGrid, LonGrid, MagnitudeCompleteness, 'LineStyle', 'none','FaceColor', 'interp')
view(90, -90)
daspect([1 LonToLatRatio 1])
h = colorbar;
set(h, 'Position', [.92 .3 .02 .4])
ylabel(h, 'Magnitude of Completeness')
daspect([1 1 1])

if ~isempty(PlotLimits)
    xlim([PlotLimits(1) PlotLimits(3)])
    ylim([PlotLimits(2) PlotLimits(4)])
else
    ylim([min(EventLon) max(EventLon)])
    xlim([min(EventLat) max(EventLat)])
end

zlim([.1 10])
caxis([min(MagnitudeCompleteness(MagnitudeCompleteness>.1)) max(max(MagnitudeCompleteness))])
disp(' ')
disp(['Circle radius: ' num2str(MaxDist)])
disp(['Grid step: ' num2str(GridStep)])
disp(['Minimum number of events: ' num2str(MinNumberOfEvents)])
disp(' ')
hold on
plot3(EventLat, EventLon, ones(length(EventLon), 1) - 3, 'o', 'Color', Grey, 'MarkerSize', 2)

title('Magnitude of Completeness (MAXC technique)', 'FontSize', 13, 'FontWeight', 'bold')
ylabel('Longitude')
xlabel('Latitude')

print(gcf,'CurrentFigures/QC_03_03_McSpaceVariation','-dpng', '-r300')








