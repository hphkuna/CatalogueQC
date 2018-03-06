%% --------
% QC_01_07_TimeVsMagCombPlot(MeasCat)

% Time Vs Magnitude plot - Comb Plot

% Input: MeasCat - Measured catalogue


function QC_01_07_TimeVsMagCombPlot(MeasCat)

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

%% display plot

% color definition
FirstColor = [0 .47 .95];
SecondColor = [.95 .47 0];
ThirdColor = [.33 .66 0];
Grey = [.7 .7 .7];

% plot
figure('name', 'Comb Plot', 'Position', [100, 100, 1049, 895])
bar(OriginTime, EventMag, 'BaseValue', 0, 'FaceColor', FirstColor)
datetick('x', 2)
xlim([min(OriginTime) max(OriginTime)])
ylim([0 1.1*max(EventMag)])

title('Comb Plot', 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Date')
ylabel('Magnitude')

print(gcf,'CurrentFigures/QC_01_07_TimeVsMagCombPlot','-dpng', '-r300')

