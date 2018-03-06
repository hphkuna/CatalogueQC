%% --------
% QC_03_02_McTimeVariation(MeasCat)

% Reveals temporal variation of Magnitude of Completeness

% Using Maximum Curvature Technique (MAXC)
% Shifting 500 samples window in time

% Input: MeasCat - Measured catalogue


function QC_03_02_McTimeVariation(MeasCat)

close all

if ~isnumeric(MeasCat)
    error('Load catagoue in the correct format (see readme for description)')
end

%% definition of the smoothing filter

disp('The Window Length or the length of the Low-pass filter length might be changed')
WindowLength = 70;
FilterLength = 200; % must be even number

%% definition of variables

% definition of event origin time, latitude, longitude, depth and magnitude
OriginTime = datenum(MeasCat(:,1), MeasCat(:,2), MeasCat(:,3), MeasCat(:,4), MeasCat(:,5), MeasCat(:,6));
EventLat = MeasCat(:, 7);
EventLon = MeasCat(:, 8);
EventDepth = MeasCat(:, 9);
EventMag = MeasCat(:, 10);

% freeing up the memory
clear MeasCat

%% Mc computation

% definition of the shifting window
NumberOfSteps = length(OriginTime)-WindowLength;

if length(OriginTime) < WindowLength
    error('Number of events in the catalogue is lower than minimum number of events needed for the analysis')
end

% definition of MagBins
MinMag = (floor(min(EventMag)*10))/10;
MaxMag = (ceil(max(EventMag)*10))/10;
MagBins = MinMag:.1:MaxMag;

% repmat the magnitude vector
EventMagMatrix = repmat(EventMag, 1, NumberOfSteps);

% create a matrix with indexes which selects magnitudes from EventMagMatrix
SelectionBasic = repmat((1:1:WindowLength)', 1, NumberOfSteps);
Increment1 = repmat((0:1:NumberOfSteps-1), WindowLength, 1);
Increment2 = repmat((0:length(EventMag):(NumberOfSteps-1)*length(EventMag)), WindowLength, 1);
SelectionMatrix = SelectionBasic + Increment1 + Increment2;

% selected magnitudes
EventMagMatrixSelect = EventMagMatrix(SelectionMatrix);

% searching for the highest count of events in each column
NumberOfEvents = hist(EventMagMatrixSelect, MagBins, 2);
[~, MaxIndex] = max(NumberOfEvents);
MagBinsMatrix = repmat(MagBins', 1, NumberOfSteps);
MaxIndex = MaxIndex + (0:length(MagBins):(NumberOfSteps-1)*length(MagBins));
MagnitudeOfCompleteness = MagBinsMatrix(MaxIndex);

%% plot the rusults

% color definition
FirstColor = [0 .47 .95];
SecondColor = [.95 .47 0];
ThirdColor = [.33 .66 0];
Grey = [.7 .7 .7];

% lowpass filtering
ConvFilter = ones(FilterLength, 1)/FilterLength;
MagnitudeOfCompletenessFiltered = conv(MagnitudeOfCompleteness, ConvFilter, 'valid');

% plot
figure('name', 'Time Variation of Magnitude of Completeness', 'Position', [100, 100, 1049, 895])
plot(OriginTime(1:end-WindowLength), MagnitudeOfCompleteness, 'Color', Grey)
datetick('x', 2)
xlim([min(OriginTime) max(OriginTime(end-WindowLength))])
ylim([min(MagnitudeOfCompleteness)-.8 max(MagnitudeOfCompleteness)+.8])
hold on
plot(OriginTime(FilterLength/2:end-WindowLength-FilterLength/2), MagnitudeOfCompletenessFiltered, 'Color', FirstColor, 'LineWidth', 2)

title('Time Variation of Magnitude of Completeness (Lowpass filtered)', 'FontSize', 13, 'FontWeight', 'bold')
xlabel('Date')
ylabel('Magnitude of Completeness (MAXC)')

print(gcf,'CurrentFigures/QC_03_02_McTimeVariation','-dpng', '-r300')


end









