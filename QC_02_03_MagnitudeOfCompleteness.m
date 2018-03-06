%% --------
% QC_02_03_MagnitudeOfCompleteness(MeasCat)

% Computes and displays the magnitude of completeness of the measured catalogue

% Based on methods described in Mignan, A., J. Woessner (2012), Estimating the magnitude of completeness for
% earthquake catalogs, Community Online Resource for Statistical Seismicity Analysis,
% doi:10.5078/corssa-00180805. Available at http://www.corssa.org.


function QC_02_03_MagnitudeOfCompleteness(MeasCat, varargin)

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

%% definition of general parameters of catalogue MeasCat and magnitude histogram

% limits of magnitude histogram (magnitude range of the catalogue)
minMag = min(EventMag);
maxMag = max(EventMag);

% limits rounded to 1/10th of magnitude point 
minMag = floor(minMag*10)/10;
maxMag = round(maxMag*10)/10;

% definition of magnitude bins for the histogram
MagBins = minMag:.1:maxMag;

% incremental magnitude histogram
[NumberOfEvents, MagBins] = hist(EventMag, MagBins);
NumberOfEventsLog = log10(NumberOfEvents);

% cummulative magnitude histogram
CummulNumberOfEvents = fliplr(cumsum(fliplr(NumberOfEvents)));
CummulNumberOfEventsLog = log10(CummulNumberOfEvents);

%% 01 Maximum Curvature (MAXC) technique

% Wyss et al. 1999; Wiemer and Wyss 2000

% Pros: fast and straightforward, statistically robust, non-parametric
% Cons: underestimates Magnitude Of Completeness (Mc) in bulk data

MagnOfCompl_MAXC = MagBins(find(NumberOfEvents == max(NumberOfEvents)));

if numel(MagnOfCompl_MAXC)>1
    MagnOfCompl_MAXC = MagnOfCompl_MAXC(end);
end

%% 02 Goodness-of-Fit Test (GFT)

% Wiemer and Wyss 2000

% Pros: G-R deviation definition
% Cons: 90% conf. not always reached - may underestimate Mc

% computes a difference between the observed and synthetic model based on the given a and b values
% Algorithm:
% step1. Choose an initial Cutoff magnitude
% step2. Obtain a and b-value for the given Cutoff mag. by the maximum likelihood
% step3. Simulate dataset with given b-value (as a straight line from logN = a - bm) 
% step4. Compute R parameter
% step5. Mc is the lowest Cutoff magnitude with R <= 5 (for 95%) or 10 (90%)

% define Cutoff magnitude
CutoffMagnitude = minMag:.1:maxMag;
NumberOfBins = numel(CutoffMagnitude);

% definition of output matrixes
Residual = zeros(NumberOfBins, 1);

for binNo = 1:NumberOfBins
    
    % step1: constraining magnitudes larger than Cutoff magnitude
    CutoffMagnitudeForLoop = CutoffMagnitude(binNo);
    EventMagForLoopUpper = EventMag(EventMag>=CutoffMagnitude(binNo));
    EventNoForLoopUpper = numel(EventMagForLoopUpper);
    
    % step2
    bValue = (1/(sum(EventMagForLoopUpper/EventNoForLoopUpper) - CutoffMagnitude(binNo))) * log10(exp(1));
    aValue = log10(EventNoForLoopUpper) + bValue*CutoffMagnitudeForLoop;
    
    % step3
    SynthNumberOfEvents_log = aValue - bValue*CutoffMagnitude;
    
    % step4
    Residual(binNo) = ((sum(abs(CummulNumberOfEventsLog(binNo:end)-SynthNumberOfEvents_log(binNo:end)))/sum(CummulNumberOfEventsLog(binNo:end)))*100);
   
end

% step5: GFT Mc is the first cutoff magnitude where Residual <= 10;
ResidUnder10 = find(Residual<=10);

if isempty(ResidUnder10)
    warning('Residuals dont decrease under 10 percent')
    ResidUnder10 = find(Residual == min(Residual));
end

MagnOfCompl_GFT = CutoffMagnitude(ResidUnder10(1));

% color definition
FirstColor = [0 .47 .95];
SecondColor = [.95 .47 0];
ThirdColor = [.33 .66 0];
Grey = [.7 .7 .7];

% plot the GFT
figure('name', '02 Goodness-of-Fit Test (GFT)', 'Position', [100, 100, 1049, 895])
plot(CutoffMagnitude, Residual, 'o', 'Color', FirstColor, 'MarkerFaceColor', FirstColor)
xlim([minMag-.1 maxMag+.1])
ylim([0 40])
hold on
plot([minMag-.1 maxMag+.1], [10 10], '--', 'Color', Grey, 'LineWidth', 1)
plot([minMag-.1 maxMag+.1], [5 5], '--', 'Color', Grey, 'LineWidth', 1)
plot([MagnOfCompl_GFT MagnOfCompl_GFT], [0 40], '--', 'Color', SecondColor, 'LineWidth', 1)
title('02 Goodness-of-Fit Test (GFT)', 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Cutoff Magnitude')
ylabel('Residual in %')

print(gcf,'CurrentFigures/QC_02_03_MagnitudeOfCompleteness_01GFT','-dpng', '-r300')


%% 03 Mc by b-value stability (MBS)

% Cao and Gao 2002

% Pros: based on b-value stability
% Cons: may overestimate Mc, relatively high uncertainty

% if Mco (Cutoff magnitude) < Mc (Magnitude of completeness) b-value is
% underestimated, if Mco > Mc, b value forms plateau

% Algorithm:
% kind of similar as the previous technique
% step1. Choose an initial Cutoff magnitude
% step2. Calculate b using maximum likelihood method
% step3. Calculate bValueConf - confidence interval of MLE b-value
% step4. Calculate bValueAverage - average b-value over interval .5
% step5. Calculate bValueDelObs - difference between abs(bValueAverage -
%        MLE)
% step6. Mc is the first Mco at which bValueDelObs < confidence interval bValueConf

% define basic variables
bValueStep = .1;
bValueAverageOver = .5;

% define Cutoff magnitude
CutoffMagnitude = minMag:.1:maxMag;
NumberOfBins = numel(CutoffMagnitude);

% definition of output matrixes
bValue = zeros(NumberOfBins, 1);
bValueConf = zeros(NumberOfBins-bValueAverageOver/bValueStep, 1);
bValueDelObs = zeros(NumberOfBins-bValueAverageOver/bValueStep, 1);
bValueAverage = zeros(NumberOfBins-bValueAverageOver/bValueStep, 1);


for binNo = 1:NumberOfBins
    
    % step1: constraining magnitudes larger than Cutoff magnitude
    CutoffMagnitudeForLoop = CutoffMagnitude(binNo);
    EventMagForLoopUpper = EventMag(EventMag>=CutoffMagnitudeForLoop);
    EventNoForLoopUpper = numel(EventMagForLoopUpper);
    
    % step2
    bValue(binNo) = (1/(sum(EventMagForLoopUpper/EventNoForLoopUpper) - CutoffMagnitudeForLoop)) * log10(exp(1));
    
end
    

for binNo = 1:NumberOfBins-bValueAverageOver/bValueStep
    
    % step1 (redefinition for the loop): constraining magnitudes larger than Cutoff magnitude
    CutoffMagnitudeForLoop = CutoffMagnitude(binNo);
    EventMagForLoopUpper = EventMag(EventMag>=CutoffMagnitudeForLoop);
    EventNoForLoopUpper = numel(EventMagForLoopUpper);
    
    % step3
    AverageEventMagnitude = mean(EventMagForLoopUpper);
    bValueConf(binNo) = 2.3*bValue(binNo)^2 * sqrt(sum((EventMagForLoopUpper-AverageEventMagnitude).^2)/(EventNoForLoopUpper*(EventNoForLoopUpper-1)));
    if isnan(bValueConf(binNo)) == 1
        bValueConf(binNo) = bValueConf(binNo-1);
    end
        
    % step4
    bValueAverage(binNo) = sum(bValue(binNo:binNo+5)) * bValueStep/(bValueAverageOver+bValueStep);
    
    % step5
    bValueDelObs(binNo) = abs(bValueAverage(binNo) - bValue(binNo)); 
    
end

% step6.
ObsLowerThanConfidence = find(bValueDelObs < bValueConf);
MagnOfCompl_MBS = CutoffMagnitude(ObsLowerThanConfidence(1));

if isempty(MagnOfCompl_MBS)
    disp('MBS Magnitude of Completeness cannot be determined')
end
    
% plot the MBS
figure('name', '03 Mc by b-value stability (MBS)', 'Position', [100, 100, 1049, 895])
patch([CutoffMagnitude(1:NumberOfBins-bValueAverageOver/bValueStep) ...
    fliplr(CutoffMagnitude(1:NumberOfBins-bValueAverageOver/bValueStep))], ...
    [bValue(1:NumberOfBins-bValueAverageOver/bValueStep) + bValueConf; ...
    flipud(bValue(1:NumberOfBins-bValueAverageOver/bValueStep) - bValueConf)], Grey, ... 
    'EdgeColor', 'none')
xlim([minMag-.1 maxMag+.1])
ylim([0 2])
hold on
plot(CutoffMagnitude, bValue, ... 
    '-o', 'Color', FirstColor, 'MarkerFaceColor', FirstColor)
plot(CutoffMagnitude(1:NumberOfBins-bValueAverageOver/bValueStep), bValueAverage, ...
    '-o','MarkerFaceColor', SecondColor, 'Color', SecondColor)
plot([MagnOfCompl_MBS MagnOfCompl_MBS], [0 2], '--','Color', ThirdColor)

legend('MLE error limit', 'MLE b-value', 'Average b-value', 'Location','northwest')
xlabel('Cutoff Magnitude')
ylabel('b-value')
title('03 Mc by b-value stability (MBS)', 'FontSize', 16, 'FontWeight', 'bold')

print(gcf,'CurrentFigures/QC_02_03_MagnitudeOfCompleteness_02MBS','-dpng', '-r300')


%% 04 Mc from the Entire Magnitude Range (OK1993, EMR)

% Woessner and Wiemer 2005

% Pros: complete frequency-magnitude distribution (FMD) model
% Cons: assumption below Mc, 4 parameters to fit

% The method is based on interpolation of probability distribution function
% in the whole interval of magnitudes
% The probability distribution function used for testing is a product of
% cummulative normal distribution attenuated by exponential decrease
% This testing function is interpolated in incremental FMD
% Algorithm:
% step1. Define the probability distribution function
% step2. Define grid for parameter search
% step3. Compute all the probability distribution functions in grid
% step4. Compute errors of probab. distributions and find the one with the
%        least square error
% step5. Mc is defined implicitly as mu + n*sigma, where mu and sigma are
%        parameters of cumulative normal distribution and n depends on how save we
%        want to be about the Mc -> n = 1: 68%; n = 2: 95%; n = 3: 99%  


% disribution definition
% step1. product of exponential and cummulative normal distribution
ProbabDensityFct = @(MagBins,b,mu2,sigma2) ...
                         (exp(-b*log(10)*MagBins) .* normcdf(MagBins,mu2,sigma2));

% normalized number of events so the histogram fits the probability distribution
NumberOfEventsNorm = NumberOfEvents/sum(NumberOfEvents);

% step2. definition of grid search range
bValueRange = .1:.05:2;
muRange = 0:.05:5;
sigmaRange = .1:.01:2;

% definition af the grid
ParameterGrid = combvec(bValueRange, muRange, sigmaRange)';

% definition of the output matrix where all the probability density functions are going to be stored
ProbabilityFctDiff = zeros(size(ParameterGrid,1), numel(MagBins));

% step3. computing probability level for all grid points (loop over magnitude bins)
for GridPoint = 1:numel(MagBins)
    ProbabilityFctDiff(:, GridPoint) = ProbabDensityFct(MagBins(GridPoint), ParameterGrid(:,1), ParameterGrid(:,2), ParameterGrid(:,3));
end

% normalizing the probability density functions so the integral of each
% function is = 1
ProbabFunctionsNorm = ProbabilityFctDiff.*repmat(1./sum(ProbabilityFctDiff,2), 1, size(ProbabilityFctDiff,2));

% step 4. computing errors of the probability functions to the real distribution
ProbabilityFctDiff = (ProbabFunctionsNorm - repmat(NumberOfEventsNorm, size(ParameterGrid,1), 1)).^2;
ProbabilityFctDiff = sqrt(sum(ProbabilityFctDiff,2));

% drumrool please... the best probability distribution!
TheRightOne = find(ProbabilityFctDiff == min(ProbabilityFctDiff));

disp(' ')
disp('04 EMR Probability Distribution Function parameters')
disp(' ')

if numel(TheRightOne)>1
    disp('Warning: The solution of the EMR test is not unique')
    TheRightOne = TheRightOne(1);
end

% and final parameters!
bValueEMR = ParameterGrid(TheRightOne,1);
muValueEMR = ParameterGrid(TheRightOne,2);
sigmaValueEMR = ParameterGrid(TheRightOne,3);

if isempty(bValueEMR)
    bValueEMR = NaN;
end
if isempty(muValueEMR)
    muValueEMR = NaN;
end
if isempty(sigmaValueEMR)
    sigmaValueEMR = NaN;
end

MagnOfCompl_EMR68 = muValueEMR + sigmaValueEMR;
MagnOfCompl_EMR95 = muValueEMR + 2*sigmaValueEMR;
MagnOfCompl_EMR99 = muValueEMR + 3*sigmaValueEMR;

% display the interpolated values
disp(' ')
display(['b-Value: ' num2str(bValueEMR)])
display(['mu: ' num2str(muValueEMR)])
display(['sigma: ' num2str(sigmaValueEMR)])
disp(' ')

% plot the results
figure('name', '04 Mc from the Entire Magnitude Range (EMR)', 'Position', [100, 100, 1049, 895])
plot(MagBins, log10(ProbabFunctionsNorm(TheRightOne,:)), 'Color', SecondColor)
hold on
plot(MagBins, log10(NumberOfEventsNorm), 'o', 'Color', FirstColor, 'MarkerFaceColor', FirstColor)
plot(MagBins, log10(CummulNumberOfEvents/CummulNumberOfEvents(1)), 'o', 'Color', Grey, 'MarkerFaceColor', Grey)
plot([MagnOfCompl_EMR68 MagnOfCompl_EMR68], [-5 1], '--', 'Color', ThirdColor)
plot([MagnOfCompl_EMR95 MagnOfCompl_EMR95], [-5 1], '--', 'Color', ThirdColor)
plot([MagnOfCompl_EMR99 MagnOfCompl_EMR99], [-5 1], '--', 'Color', ThirdColor)
xlim([minMag-.1 maxMag+.1])

xlabel('Magnitude')
ylabel('Normalized Number of Events')
title('04 Mc from the Entire Magnitude Range (EMR)', 'FontSize', 16, 'FontWeight', 'bold')
legend('Best fit probability distribution function', 'Incremental FMD', 'Comulative FMD','Left to right: 68%, 95%, 99% Mc')

print(gcf,'CurrentFigures/QC_02_03_MagnitudeOfCompleteness_03EMR','-dpng', '-r300')


%% 05 Median-Based Analysis of the Segment Slope (MBASS)

% Amorese 2007

% Pros: non-parametric
% Cons: main discontinuity may not be Mc, relatively high uncertainty

%% Print results

disp(' ')
disp('---------------')
disp('02_03 Magnitude of Completeness - Output:')
disp(' ')
disp(['Maximum Curvature (MAXC) technique: ' num2str(MagnOfCompl_MAXC)])
disp(['Goodness-of-Fit Test (GFT): ' num2str(MagnOfCompl_GFT)])
disp(['Mc by b-value stability (MBS): ' num2str(MagnOfCompl_MBS)])
disp(['Mc from the Entire Magnitude Range (EMR) - 68%: ' num2str(MagnOfCompl_EMR68)])
disp(['Mc from the Entire Magnitude Range (EMR) - 95%: ' num2str(MagnOfCompl_EMR95)])
disp(['Mc from the Entire Magnitude Range (EMR) - 99%: ' num2str(MagnOfCompl_EMR99)])
disp('---------------')
disp(' ')

%% Plot results

figure('name', 'Magnitude of Completeness computed using different methods', 'Position', [100, 100, 1049, 895])

plot(MagBins, NumberOfEventsLog, 'o', 'Color', FirstColor, 'MarkerFaceColor', FirstColor)
hold on
plot(MagBins, CummulNumberOfEventsLog, 'o', 'Color', SecondColor, 'MarkerFaceColor', SecondColor)
plot([MagnOfCompl_MAXC MagnOfCompl_MAXC], [0 5], '--', 'Color', [.8 .4 0])
plot([MagnOfCompl_GFT MagnOfCompl_GFT], [0 5], '--', 'Color', [0 .4 .8])
plot([MagnOfCompl_MBS MagnOfCompl_MBS], [0 5], '--', 'Color', [.8 0 .4])
plot([MagnOfCompl_EMR68 MagnOfCompl_EMR68], [0 5], '--', 'Color', [.3 .6 0])
plot([MagnOfCompl_EMR95 MagnOfCompl_EMR95], [0 5], '--', 'Color', [.3 .6 0])

if ~isempty(PlotLimits)
    xlim([PlotLimits(1) PlotLimits(3)])
    ylim([PlotLimits(2) PlotLimits(4)])
else
    xlim([minMag-.1 maxMag+.1])
    ylim([0 1.1*max(CummulNumberOfEventsLog)])
end

title('Magnitude of Completeness computed using different methods', 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Magnitude')
ylabel('Number of Events')
legend('Incremental FMD', 'Cumulative FMD', 'MAXC', 'GFT', 'MBS', 'EMR68','EMR95')

print(gcf,'CurrentFigures/QC_02_03_MagnitudeOfCompleteness_04all','-dpng', '-r300')

















