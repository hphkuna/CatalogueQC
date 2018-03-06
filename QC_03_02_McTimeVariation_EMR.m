%% --------
% QC_03_02_McTimeVariation(MeasCat)

% Reveals temporal variation of Magnitude of Completeness

% Using Maximum Curvature Technique (MAXC)
% Shifting 500 samples window in time

% Input: MeasCat - Measured catalogue


function QC_03_02_McTimeVariation_EMR(MeasCat)

close all

if ~isnumeric(MeasCat)
    error('Load catagoue in the correct format (see readme for description)')
end

%% Set limits

OriginTime = datenum(MeasCat(:,1), MeasCat(:,2), MeasCat(:,3), MeasCat(:,4), MeasCat(:,5), MeasCat(:,6));

Step = 5;
WindowLength = 30;
MinOrigTime = min(OriginTime);
MaxOrigTime = max(OriginTime);

OrigTimeMins = MinOrigTime : Step : MaxOrigTime-WindowLength;
OrigTimeMaxs = OrigTimeMins + WindowLength;

NumberOfSteps = numel(OrigTimeMins)

%% Main loop

MagnCompl = zeros(NumberOfSteps, 1);

for n = 1:NumberOfSteps
    
    disp(n)
    
    MeasCatSub = MeasCat(OriginTime>OrigTimeMins(n) & OriginTime<OrigTimeMaxs(n), :);
    
    MagnCompl(n) = EMR(MeasCatSub);
    
end

%% plot the rusults

% color definition
FirstColor = [0 .47 .95];
SecondColor = [.95 .47 0];
ThirdColor = [.33 .66 0];
Grey = [.7 .7 .7];
 
% plot
figure('name', 'Time Variation of Magnitude of Completeness', 'Position', [100, 100, 1049, 895])
plot(OrigTimeMins, MagnCompl, 'Color', SecondColor, 'LineWidth', 2)
datetick('x', 2)
xlim([min(OriginTime) max(OriginTime(end-WindowLength))])
ylim([min(MagnCompl)-.8 max(MagnCompl)+.8])


title('Time Variation of Magnitude of Completeness (Lowpass filtered)', 'FontSize', 13, 'FontWeight', 'bold')
xlabel('Date')
ylabel('Magnitude of Completeness (EMR)')

print(gcf,'CurrentFigures/QC_03_02_McTimeVariation_EMR','-dpng', '-r300')


end

function MC_EMR68 = EMR(MeasCat)


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

MC_EMR68 = muValueEMR + sigmaValueEMR;

end








