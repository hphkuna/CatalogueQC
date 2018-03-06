%% --------
% QC_03_01_ArtificialEqTest(MeasCat)

% The procedure tests a presence of man-induced earthquakes (explosions etc)

% Input: MeasCat - Measured catalogue, Mc - Magnitude of Completesess


function QC_03_01_ArtificialEqTest(MeasCat, Mc)

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

EventOriginHour = MeasCat(:,4);

% freeing up the memory
clear MeasCat

%% display plots

% color definition
FirstColor = [0 .47 .95];
SecondColor = [.95 .47 0];
ThirdColor = [.33 .66 0];
Grey = [.7 .7 .7];

% Hour of the days test
HourBins = 0:1:23;

figure('name', 'Hour of the day test', 'Position', [100, 100, 1300, 500])
subplot(1,2,1)
hist(EventOriginHour, HourBins)
h = findobj(gca,'Type','patch');
set(h, 'FaceColor', FirstColor)
set(h, 'EdgeColor', 'w');

xlim([-.5 23.5])

title('Hour of the day test - All magnitudes', 'FontSize', 11, 'FontWeight', 'bold')
xlabel('Hour of the day')
ylabel('Number of events')

subplot(1,2,2)
hist(EventOriginHour(EventMag>=Mc), HourBins)
h = findobj(gca,'Type','patch');
set(h, 'FaceColor', FirstColor)
set(h, 'EdgeColor', 'w');

xlim([-.5 23.5])

title('Hour of the day test - magnitudes greater than Mc', 'FontSize', 11, 'FontWeight', 'bold')
xlabel('Hour of the day')
ylabel('Number of events')

print(gcf,'CurrentFigures/QC_03_01_ArtificialEqTest_Hour','-dpng', '-r300')

% Day of the week test
DayNumber = weekday(OriginTime)+1;
DayBins = 2:1:8;

figure('name', 'Day of the week test', 'Position', [100, 100, 1300, 500])
subplot(1,2,1)
hist(DayNumber, DayBins)
h = findobj(gca,'Type','patch');
set(h, 'FaceColor', FirstColor)
set(h, 'EdgeColor', 'w');

datetick('x', 'ddd')
xlim([1.5 8.5])

title('Day of the week test - All magnitudes', 'FontSize', 11, 'FontWeight', 'bold')
xlabel('Day of the week')
ylabel('Number of events')

subplot(1,2,2)
hist(DayNumber(EventMag>=Mc), DayBins)
h = findobj(gca,'Type','patch');
set(h, 'FaceColor', FirstColor)
set(h, 'EdgeColor', 'w');

datetick('x', 'ddd')
xlim([1.5 8.5])

title('Day of the week test - magnitudes greater than Mc', 'FontSize', 11, 'FontWeight', 'bold')
xlabel('Day of the week')
ylabel('Number of events')

print(gcf,'CurrentFigures/QC_03_01_ArtificialEqTest_Day','-dpng', '-r300')

