%% --------
% QC_02_01_MagHist(MeasCat)

% Displays histogram and cummulative histogram of earthquake magnitudes
% Prints minimum and maximum magnitude, 10 largest events in the catalogue

function QC_02_01_MagHist(MeasCat, varargin)

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

%% histogram of events

% limits of magnitude histogram (magnitude range of the catalogue)
minMag = min(EventMag);
maxMag = max(EventMag);

% limits rounded to 1/10th of magnitude point 
minMag = floor(minMag*10)/10;
maxMag = ceil(maxMag*10)/10;

% definition of magnitude bins for the histogram
MagBins = minMag:.1:maxMag;

% incremental magnitude histogram
[NumberOfEvents, MagBins] = hist(EventMag, MagBins);
NumberOfEventsLog = log10(NumberOfEvents);

% cummulative magnitude histogram
CummulNumberOfEvents = fliplr(cumsum(fliplr(NumberOfEvents)));
CummulNumberOfEventsLog = log10(CummulNumberOfEvents);

%% display of Magnitude histogram

% color definition
FirstColor = [0 .47 .95];
SecondColor = [.95 .47 0];
ThirdColor = [.33 .66 0];
Grey = [.7 .7 .7];

% incremental histogram
figure('name', 'Frequency-Magnitude Distribution', 'Position', [100, 100, 1049, 895])
plot(MagBins, NumberOfEventsLog, 'o', 'Color', FirstColor, 'LineWidth', 1, 'MarkerFaceColor', FirstColor)
hold on

if ~isempty(PlotLimits)
    xlim([PlotLimits(1) PlotLimits(3)])
    ylim([PlotLimits(2) PlotLimits(4)])
else
    xlim([minMag-.1 maxMag+.1])
end

xlabel('Event Magnidute')
ylabel('Number of Events')
% cummulative histogram
plot(MagBins, CummulNumberOfEventsLog, 'o', 'Color', SecondColor, 'LineWidth', 1, 'MarkerFaceColor', SecondColor)

title('Frequency-magnitude distribution', 'FontSize', 16, 'FontWeight', 'bold')
legend('incremental histogram of magnitude distribution', 'cummulative histogram of magnitude distribution')

print(gcf,'CurrentFigures/QC_02_01_MagHist','-dpng', '-r300')

%% print values

disp(' ')
disp('---------------')
disp('02_01 Frequency-magnitude distribution - Output:')
disp(['Minimum magnitude: ' num2str(min(EventMag))])
disp(['Maximum magnitude: ' num2str(max(EventMag))])
disp(' ')

% 10 largest events
DecreasEventMag = sort(EventMag, 'descend');
disp('10 largest events in the catalogue:')
disp(num2str(DecreasEventMag(1:10)))
disp('---------------')
disp(' ')












