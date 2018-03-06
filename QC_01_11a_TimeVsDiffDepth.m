%% --------
% QC_01_10_TimeVsDiffDepth(MeasCat)

% Displays scatter plot of tive vs depth

% Input: MeasCat - Measured catalogue


function QC_01_11a_TimeVsDiffDepth(MeasCat, varargin)

FilterLength = 20;

% close all

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

%% diferential depths and smoothing

DiffEventDepth = diff(EventDepth);

filter = ones(FilterLength,1)/FilterLength;
DiffEventDepthFilt = conv(DiffEventDepth, filter, 'valid');
DiffEventDepthFilt2 = conv(DiffEventDepthFilt, filter, 'valid');

%% display plot

% color definition
FirstColor = [0 .47 .95];
SecondColor = [.95 .47 0];
ThirdColor = [.33 .66 0];
Grey = [.7 .7 .7];

% plot
MinTime = min(OriginTime);
MaxTime = max(OriginTime);

% figure('name', 'Time Vs. Diff Depth diagram', 'Position', [100, 100, 800, 800])
% subplot(2,1,1)
% plot(OriginTime(2:end), DiffEventDepth, 'o', 'Color', FirstColor, 'MarkerSize', 4)
% datetick('x',2)
% 
% if ~isempty(PlotLimits)
%     ylim([PlotLimits(1) PlotLimits(3)])
%     xlim([PlotLimits(2) PlotLimits(4)])
% else
%     xlim([MinTime MaxTime])
%     ylim([min(DiffEventDepth) max(DiffEventDepth)])
% end
% 
% title('Time Vs. Difference depth diagram', 'FontSize', 16, 'FontWeight', 'bold')
% xlabel('Time')
% ylabel('Depth difference between consecutive events')

subplot(3,1,3)
plot(OriginTime(FilterLength:end-FilterLength), DiffEventDepthFilt2, '-o','Color', SecondColor, 'MarkerSize', 2, 'LineWidth', 2)
datetick('x',2)

% StDev = std(DiffEventDepthFilt2);
StDev = 0.045;

hold on
MinTime = min(OriginTime);
MaxTime = max(OriginTime);
plot([MinTime MaxTime], [0 0], 'Color', FirstColor)
plot([MinTime MaxTime], [StDev StDev], ':', 'Color', FirstColor)
plot([MinTime MaxTime], [-StDev -StDev], ':', 'Color', FirstColor)
plot([MinTime MaxTime], [2*StDev 2*StDev], '--', 'Color', FirstColor)
plot([MinTime MaxTime], [-2*StDev -2*StDev], '--', 'Color', FirstColor)

title('Time Vs. Difference depth diagram - Low-pass filtered', 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Time')
ylabel('Depth difference between consecutive events')

if ~isempty(PlotLimits)
    ylim([PlotLimits(1) PlotLimits(3)])
    xlim([PlotLimits(2) PlotLimits(4)])
else
    xlim([MinTime MaxTime])
    ylim([min(DiffEventDepthFilt2) -min(DiffEventDepthFilt2)])
end


print(gcf,'CurrentFigures/QC_01_11_TimeVsDiffDepth','-dpng', '-r300')













