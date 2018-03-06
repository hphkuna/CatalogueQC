%% --------
% QC_01_08_CumulativeEventCount(MeasCat)

% Cumulative Event Count

% Input: MeasCat - Measured catalogue


function NumberOfEventsCumul = QC_01_08_CumulativeEventCount(MeasCat)

% close all

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

TimeBins = floor(min(OriginTime)):1:ceil(max(OriginTime));
NumberOfEvents = hist(OriginTime, TimeBins);

NumberOfEventsCumul = cumsum(NumberOfEvents);

%% display plot

% color definition
FirstColor = [0 .47 .95];
SecondColor = [.95 .47 0];
ThirdColor = [.33 .66 0];
Grey = [.7 .7 .7];

% plot
% figure('name', 'Cumulative Event Count in Time', 'Position', [100, 100, 1049, 895])
plot(TimeBins, NumberOfEventsCumul/sum(NumberOfEventsCumul), 'Color', FirstColor)
datetick('x', 2)
xlim([min(TimeBins) max(TimeBins)])

title('Cumulative Event Count in Time', 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Date')
ylabel('Number Of Events')

print(gcf,'CurrentFigures/QC_01_08_CumulativeEventCount','-dpng', '-r300')

