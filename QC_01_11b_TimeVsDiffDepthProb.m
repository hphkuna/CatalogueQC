%% --------
% QC_01_11b_TimeVsDiffDepthProb(MeasCat)

% Displays scatter plot of time vs depth

% Input: MeasCat - Measured catalogue


function QC_01_11b_TimeVsDiffDepthProb(MeasCat, varargin)

WindowLength = 20;
FilterLength = 20;

close all

if ~isnumeric(MeasCat)
    error('Load catagoue in the correct format (see readme for description)')
end

PlotLimits = cell2mat(varargin);

% MeasCat = MeasCat(MeasCat(:,9)>10,:);

%% definition of variables

% definition of event origin time, latitude, longitude, depth and magnitude
OriginTime = datenum(MeasCat(:,1), MeasCat(:,2), MeasCat(:,3), MeasCat(:,4), MeasCat(:,5), MeasCat(:,6));
EventLat = MeasCat(:, 7);
EventLon = MeasCat(:, 8);
EventDepth = MeasCat(:, 9);
EventMag = MeasCat(:, 10);

% freeing up the memory
clear MeasCat

%% main loop

NumOfEvents = numel(OriginTime);
NumOfSteps = NumOfEvents - WindowLength;

Percentile = zeros(NumOfSteps,1);

for n = 1:NumOfSteps
    
    DepthForLoop = EventDepth(n) - EventDepth(n+1:n+WindowLength);
    Percentile(n) = sum(DepthForLoop<0)/WindowLength;
        
end

%% diferential depths and smoothing

filter = ones(FilterLength,1)/FilterLength;
PercentileFilt = conv(Percentile, filter, 'valid');
PercentileFilt2 = conv(PercentileFilt, filter, 'valid');

%% display plot

% color definition
FirstColor = [0 .47 .95];
SecondColor = [.95 .47 0];
ThirdColor = [.33 .66 0];
Grey = [.7 .7 .7];

% plot
MinTime = min(OriginTime);
MaxTime = max(OriginTime(1:NumOfSteps));

figure('name', 'Time Vs. Percentile diagram', 'Position', [100, 100, 800, 800])
subplot(2,1,1)
plot(OriginTime(1:NumOfSteps), Percentile, 'o', 'Color', FirstColor, 'MarkerSize', 4)
datetick('x',2)

if ~isempty(PlotLimits)
    ylim([PlotLimits(1) PlotLimits(3)])
    xlim([PlotLimits(2) PlotLimits(4)])
else
    xlim([MinTime MaxTime])
    ylim([min(Percentile) max(Percentile)])
end

title('Time Vs. Percentile', 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Time')
ylabel('Depth difference between consecutive events')

subplot(2,1,2)
plot(OriginTime(FilterLength:NumOfSteps-FilterLength+1), PercentileFilt2, '-o','Color', SecondColor, 'MarkerSize', 2, 'LineWidth', 2)
datetick('x',2)

% StDev = std(PercentileFilt2);
StDev = .0516;

hold on
plot([MinTime MaxTime], [.5 .5], 'Color', FirstColor)
plot([MinTime MaxTime], [.5+StDev .5+StDev], ':', 'Color', FirstColor)
plot([MinTime MaxTime], [.5-StDev .5-StDev], ':', 'Color', FirstColor)
plot([MinTime MaxTime], [.5+2*StDev .5+2*StDev], '--', 'Color', FirstColor)
plot([MinTime MaxTime], [.5-2*StDev .5-2*StDev], '--', 'Color', FirstColor)

title('Time Vs. Percentile diagram - Low-pass filtered', 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Time')
ylabel('Depth difference between consecutive events')

if ~isempty(PlotLimits)
    ylim([PlotLimits(1) PlotLimits(3)])
    xlim([PlotLimits(2) PlotLimits(4)])
else
    xlim([MinTime MaxTime])
    ylim([min(PercentileFilt2) max(PercentileFilt2)])
end


print(gcf,'CurrentFigures/QC_01_11b_TimeVsPercentile','-dpng', '-r300')


