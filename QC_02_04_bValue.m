%% --------
% QC_02_04_bValue(MeasCat, Mc)

% Computes and displays the b-value of dataset

% The b-value estimated according to Aki K (1965). Maximum Likelihood Estimate of b in the Formula logN = a - bM and its confidence limits.
% Bulletin of the Earthquake Research Institute, 43, 237-239.

% Input: MeasCat - Measured catalogue, Mc - Magnitude of Completeness


function QC_02_04_bValue(MeasCat, MagnitudeCompleteness, varargin)

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

% limits of magnitude for Cutoff magnitude values
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

%% maximum likelihood estimation of b value for different Cutoff magnitude

if MagnitudeCompleteness < minMag || MagnitudeCompleteness > maxMag
    disp(['Minimum catalogue magnitude: ' num2str(minMag)])
    disp(['Maximum catalogue magnitude: ' num2str(maxMag)])
    error('Error: Mc out of bounds')
end
  
% constraining magnitudes larger than Cutoff magnitude
EventMagForMc = EventMag(EventMag>=MagnitudeCompleteness);
EventNoForMc = numel(EventMagForMc);
MagBinsBValue = MagnitudeCompleteness:.1:maxMag;

% mle estimator
myGutPDF = @(mags, bValue) xGutenbergPDF(mags, bValue);

[bValue, bValueErr95]  = mle(EventMagForMc, 'pdf', myGutPDF, 'start', .8, 'lowerbound', 0, 'upperbound', 100);
bValueErr95Int = bValueErr95(2) - bValue;

aValue = log10(max(myGutPDF(MagBinsBValue, bValue)) * EventNoForMc * .1);
aValueErr95 = zeros(2,1);
aValueErr95(1) = log10(max(myGutPDF(MagBinsBValue, bValueErr95(1))) * EventNoForMc * .1);
aValueErr95(2) = log10(max(myGutPDF(MagBinsBValue, bValueErr95(2))) * EventNoForMc * .1);
aValueErr95Int = aValueErr95(2) - aValue;

% regressed line
RegressedLine = aValue - bValue*(MagBinsBValue-MagnitudeCompleteness);

%% least squares estimation from the cummulative event count

MagBins4LS = (MagBins(MagBins>=MagnitudeCompleteness))';
CummulNumberLog4LS = (CummulNumberOfEventsLog(MagBins>=MagnitudeCompleteness))';

G = [MagBins4LS ones(numel(MagBins4LS),1)];
m = G\CummulNumberLog4LS;

CummulativeRegressedLine = m(2) + m(1)*MagBins4LS;

StanDeviation = sum((CummulNumberLog4LS - G*m).^2)/numel(CummulativeRegressedLine);
Covar = StanDeviation * inv(G'*G);
bValueCum95Int = 1.96 * sqrt(Covar(4));

%% plotting the b-value with respect to Cutoff magnitude used for calculation

% color definition
FirstColor = [0 .47 .95];
SecondColor = [.95 .47 0];
ThirdColor = [.33 .66 0];
Grey = [.7 .7 .7];

% sublot of incremental and cummulative histogram
figure('name', 'Frequency-magnitude distribution and b-value', 'Position', [100, 100, 1049, 895])
% incremental histogram
plot(MagBins, NumberOfEventsLog, 'o', 'Color', FirstColor, 'LineWidth', 1, 'MarkerFaceColor', FirstColor)
hold on
xlabel('Event Magnidute')
ylabel('Number of Events')

if ~isempty(PlotLimits)
    xlim([PlotLimits(1) PlotLimits(3)])
    ylim([PlotLimits(2) PlotLimits(4)])
else
    xlim([minMag-.1 maxMag+.1])
    ylim([0 1.1*max(CummulNumberOfEventsLog)])
end

% cummulative histogram
plot(MagBins, CummulNumberOfEventsLog, 'o', 'Color', SecondColor, 'LineWidth', 1, 'MarkerFaceColor', SecondColor)
% regressed line
plot(MagBinsBValue, RegressedLine, 'Color', ThirdColor, 'LineWidth', 2)
plot(MagBins4LS, CummulativeRegressedLine, 'Color', ThirdColor, 'LineWidth', 2)

text(MagBinsBValue(1), RegressedLine(1) - .4, num2str(bValue))
text(MagBins4LS(1), CummulativeRegressedLine(1) + .1, num2str(-m(1)))

title('Frequency-magnitude distribution', 'FontSize', 16, 'FontWeight', 'bold')
legend('incremental histogram of magnitude distribution', 'cummulative histogram of magnitude distribution')

print(gcf,'CurrentFigures/QC_02_04_bValue','-dpng', '-r300')

%% print results

disp(' ')
disp('---------------')
disp('QC_02_04 bValue - Output:')
disp(' ')
disp(['Magnitude of Completeness: ' num2str(MagnitudeCompleteness)])
disp(' ')
disp(['bValue: ' num2str(bValue)])
disp(['bValue 95% confidence: ' num2str(bValueErr95Int)])
disp(' ')
disp(['bValue Cummulative: ' num2str(-m(1))])
disp(['bValue Cummulative 95% confidence: ' num2str(bValueCum95Int)])
disp('---------------')
disp(' ')


