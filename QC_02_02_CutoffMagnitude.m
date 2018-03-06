%% --------
% QC_02_02_CutoffMagnitude(MeasCat)

% Computes and displays the b-value with respect to cutoff magnitude

% The b-value estimated according to Aki K (1965). Maximum Likelihood Estimate of b in the Formula logN = a - bM and its confidence limits.
% Bulletin of the Earthquake Research Institute, 43, 237-239.


function QC_02_02_CutoffMagnitude(MeasCat)

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

%% maximum likelihood estimation of b value for different Cutoff magnitude

% limits of magnitude for Cutoff magnitude values
minMag = min(EventMag);
maxMag = max(EventMag);

% limits rounded to 1/10th of magnitude point 
minMag = floor(minMag*10)/10;
maxMag = ceil(maxMag*10)/10;

% definition of Cutoff magnitude
CutoffMagnitude = minMag:.1:maxMag;

% number of events
NumberOfEvents = numel(EventMag);

% number of magnitude bins
NumberOfBins = numel(CutoffMagnitude);

% definition of the output matrixes
bValue = zeros(NumberOfBins, 1);
bValue_low95lim = zeros(NumberOfBins, 1);
bValue_up95lim = zeros(NumberOfBins, 1);

% D_eps values for the confidence error bars
D_eps_50 = .66;
D_eps_80 = 1.30;
D_eps_90 = 1.64;
D_eps_95 = 1.96;
D_eps_98 = 2.34;

% loop computing the b value and the upper and lower confidence limits by
% maximum likelihood method
for binNo = 1:NumberOfBins
    
    % constraining magnitudes larger than Cutoff magnitude
    EventMagForLoop = EventMag(EventMag>=CutoffMagnitude(binNo));
    EventNoForLoop = numel(EventMagForLoop);

    % calculating b values and confidence limits
    bValue(binNo) = (1/(sum(EventMagForLoop/EventNoForLoop) - CutoffMagnitude(binNo))) * log10(exp(1));
    bValue_low95lim(binNo) = ((1 - D_eps_95/sqrt(EventNoForLoop))/(sum(EventMagForLoop/EventNoForLoop) - CutoffMagnitude(binNo))) * log10(exp(1));
    bValue_up95lim(binNo) = ((1 + D_eps_95/sqrt(EventNoForLoop))/(sum(EventMagForLoop/EventNoForLoop) - CutoffMagnitude(binNo))) * log10(exp(1));
    
end

bValue_low95lim = abs(bValue_low95lim - bValue);
bValue_up95lim = abs(bValue_up95lim - bValue);

%% histogram of events

% definition of magnitude bins for the histogram
MagBins = minMag:.1:maxMag;

% incremental magnitude histogram
[NumberOfEvents, MagBins] = hist(EventMag, MagBins);
NumberOfEventsLog = log10(NumberOfEvents);

% cummulative magnitude histogram
CummulNumberOfEvents = fliplr(cumsum(fliplr(NumberOfEvents)));
CummulNumberOfEventsLog = log10(CummulNumberOfEvents);

%% plotting the b-value with respect to Cutoff magnitude used for calculation

% color definition
FirstColor = [0 .47 .95];
SecondColor = [.95 .47 0];
ThirdColor = [.33 .66 0];
Grey = [.7 .7 .7];

% sublot of incremental and cummulative histogram
figure('name', 'Frequency-magnitude distribution and b-value', 'Position', [100, 100, 1049, 895])

% subplot of b Value estimation
errorbar(CutoffMagnitude, bValue, bValue_low95lim, bValue_up95lim, 'o', 'Color', FirstColor, 'LineWidth', 1, 'MarkerFaceColor', FirstColor)
xlabel('Cutoff Magnitude')
ylabel('b-value Estimation')
xlim([minMag-.1 maxMag+.1])
ylim([-1 6])

title('b-value respective to the magnitude cutoff', 'FontSize', 16, 'FontWeight', 'bold')
legend('b-value from the incremental histogram', 'Location','northwest')

print(gcf,'CurrentFigures/QC_02_02_CutoffMagnitude','-dpng', '-r300')

    
    
    
    
    
    
    
    
    
    
    