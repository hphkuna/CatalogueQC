%% --------
% QC_02_05_TimeVsMagnitude(MeasCat)

% Displays a scatter plot of time vs magnitude of event

% Reflects variation is seismicity (eg. aftershock sequences) as well as
% changes of network configuration

% Input: MeasCat - Measured catalogue


function QC_02_05_TimeVsMagnitude(MeasCat)

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

%% scatter plot

% color definition
FirstColor = [0 .47 .95];
SecondColor = [.95 .47 0];
ThirdColor = [.33 .66 0];
Grey = [.7 .7 .7];

figure('name', 'Date vs. magnitude plot', 'Position', [100, 100, 1049, 895])
plot(OriginTime, EventMag, 'o', 'Color', FirstColor, 'MarkerSize', 2)
datetick('x',2)
xlim([min(OriginTime) max(OriginTime)])
xlabel('Date')
ylabel('Magnitude')
title('Date vs. magnitude plot', 'FontSize', 16, 'FontWeight', 'bold')

print(gcf,'CurrentFigures/QC_02_05_TimeVsMagnitude','-dpng', '-r300')


