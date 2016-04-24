% Create RateSpec from the interest rate term structure
StartDates = '24-April-2014';
EndDates = '24-April-2020';
Rates = 0.035;
Compounding = -1;
Basis = 1;

RateSpec = intenvset('ValuationDate', StartDates, 'StartDates', StartDates, ...
    'EndDates', EndDates, 'Rates', Rates, 'Compounding', ...
    Compounding, 'Basis', Basis);

% Define StockSpec with the underlying asset information
Sigma = 0.20;
AssetPrice = 100;

StockSpec = stockspec(Sigma, AssetPrice);

% Define the Asian option
Settle = '24-April-2014';
ExerciseDates = '24-April-2015';
Strike = 100;
OptSpec = 'call';

% Levy model approximation
PriceLevy = asianbylevy(RateSpec, StockSpec, OptSpec, Strike, Settle,...
                        ExerciseDates);

% Kemna-Vorst closed form model
PriceKV = asianbykv(RateSpec, StockSpec, OptSpec, Strike, Settle,...
                     ExerciseDates);

% Comparison of calculated prices for the geometric and arithmetic options
% using different closed form algorithms.
fprintf('Comparison of Asian Arithmetic and Geometric Prices:\n\n');
fprintf('Levy:         %f\n', PriceLevy);
fprintf('Kemna-Vorst:  %f\n', PriceKV);

% Simulation Parameters
NumTrials = 500;
NumPeriods = 200;

% Price the arithmetic option
PriceAMC = asianbyls(RateSpec, StockSpec, OptSpec, Strike, Settle,...
                     ExerciseDates,'NumTrials', NumTrials, ...
                     'NumPeriods', NumPeriods);

% Price the geometric option
PriceGMC = asianbyls(RateSpec, StockSpec, OptSpec, Strike, Settle,...
                     ExerciseDates,'NumTrials', NumTrials, ...
                     'NumPeriods', NumPeriods, 'AvgType', 'arithmetic');

% Use the antithetic variates method to value the options
Antithetic = true;
PriceAMCAntithetic = asianbyls(RateSpec, StockSpec, OptSpec, Strike, Settle,...
                    ExerciseDates,'NumTrials', NumTrials, 'NumPeriods',...
                    NumPeriods, 'Antithetic', Antithetic);

PriceGMCAntithetic = asianbyls(RateSpec, StockSpec, OptSpec, Strike, Settle,...
                    ExerciseDates,'NumTrials', NumTrials, 'NumPeriods',...
                    NumPeriods, 'Antithetic', Antithetic,'AvgType', 'geometric');

% Display prices
fprintf('Asian Prices using Monte Carlo Method:\n\n');
fprintf('Arithmetic Asian\n');
fprintf('Standard Monte Carlo:           %f\n', PriceAMC);
fprintf('Variate Antithetic Monte Carlo: %f\n\n', PriceAMCAntithetic);
fprintf('Geometric Asian\n');
fprintf('Standard Monte Carlo:           %f\n', PriceGMC);
fprintf('Variate Antithetic Monte Carlo: %f\n', PriceGMCAntithetic);