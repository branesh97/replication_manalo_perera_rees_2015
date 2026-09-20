clear all
clc

%% Preliminary

% Folder containing this script
folder = fileparts(mfilename('fullpath'));

% Add functions folder
addpath(fullfile(folder,'_func'))

% load data
filePath = fullfile(folder, 'data', 'data_raw.xlsx');
data_raw = readmatrix(filePath, 'Sheet', 'Raw', 'Range', 'B2:H167');
% data_raw columns:
% 1 = US real GDP
% 2 = Australian terms of trade
% 3 = Australian real GDP
% 4 = Australian trimmed mean inflation
% 5 = Australian cash rate
%     historical official-market/cash-target series pre-Jul 1998,
%     interbank overnight cash rate thereafter; quarterly average;
%     this series data construction appears important.
% 6 = Australian real TWI
% 7 = VIX


% Timeline: 1985Q1 to 2026Q2
% Convention: Q1 = year, Q2 = year+0.25, Q3 = year+0.50, Q4 = year+0.75
yearlab =(1985:1/4:2026+1/4)'; % timeline (1985Q1 to 2026Q2)
assert(size(data_raw,1) == length(yearlab)) % check

%% Extended sample: 1985Q1 to 2026Q2

sample_ext = (yearlab >= 1985) & (yearlab <= 2026 + 1/4);

raw_ext = data_raw(sample_ext,1:6);   % VIX not used in baseline model
assert(~any(isnan(raw_ext(:))))        % check
yearlab_ext = yearlab(sample_ext);

T = size(raw_ext,1);
assert(T == 166) % check

%% Model settings
% Identical to baseline replication

p = 2;                  % Manalo et al. use two lags
n_foreign = 2;          % US GDP and terms of trade
shock_er = 6;           % Real TWI is ordered last

hor = 24;               % IRF horizon: 24 quarters
hor_fevd = 40;

%% Transform variables

% Log variables
log_us_gdp  = 100*log(raw_ext(:,1));
log_tot     = 100*log(raw_ext(:,2));
log_aus_gdp = 100*log(raw_ext(:,3));
log_rtwi    = 100*log(raw_ext(:,6));

% Already in required units
trim_inf  = raw_ext(:,4);      % quarterly percentage change
cash_rate = raw_ext(:,5);      % percent, level

%% Quadratic detrending of GDP
% IMPORTANT: trends are re-estimated over the full extended sample

tt = (1:T)';
Xtrend = [ones(T,1) tt tt.^2];

% US GDP
b_us = Xtrend \ log_us_gdp;             % coef
trend_us = Xtrend*b_us;                 % fitted value
us_gdp_gap = log_us_gdp - trend_us;     % resid = cyclical dev from trend

% Australian GDP
b_aus = Xtrend \ log_aus_gdp;           % coef
trend_aus = Xtrend*b_aus;               % fitted value
aus_gdp_gap = log_aus_gdp - trend_aus;  % resid = cyclical dev from trend

%% Dataset entering aggregate VAR
y = [us_gdp_gap ...
    log_tot ...
    aus_gdp_gap ...
    trim_inf ...
    cash_rate ...
    log_rtwi];

Variables = {'US Real GDP', ...
    'Terms of Trade', ...
    'Australian Real GDP', ...
    'Trimmed Mean Inflation', ...
    'Cash Rate', ...
    'Real TWI'};

N = size(y,2);

% Deterministic controls
% same used by Manalo et al. (2015)
% - no COVID dummies yet

% (1) Inflation-targeting regime dummy
d_it = double(yearlab_ext >= 1993);

% (2) Individual GFC quarter pulse dummies
d_2008q4 = double(yearlab_ext == 2008 + 3/4);
d_2009q1 = double(yearlab_ext == 2009);
d_2009q2 = double(yearlab_ext == 2009 + 1/4);
d_2009q3 = double(yearlab_ext == 2009 + 2/4);

d = [d_it d_2008q4 d_2009q1 d_2009q2 d_2009q3];

%% Estimate extended-sample aggregate VAR

[phi,gamma,SIGMA,X,e] = olsvar_soe(y,p,d,n_foreign);

%% Identification and impulse response functions

% Cholesky identification
% e_t = B * epsilon_t
B = chol(SIGMA,'lower');

% Normalize exchange-rate shock to a 10 percent appreciation on impact
B_N = B;

B_N(:,shock_er) = 10*B(:,shock_er)./B(shock_er,shock_er);

% IRFs following normalized 10 percent appreciation
IRF_N = calculate_IRF_FEVD(phi,B_N,hor);

%% Extract responses to exchange-rate shock

ind_er = (shock_er-1)*N + (1:N);

IRF_ER = IRF_N(ind_er,:);

%% Convert quarterly inflation IRF to year-ended inflation IRF

irf_inf_q = IRF_ER(4,:);

temp = conv(irf_inf_q,ones(1,4));   % help with running sum of 4 quarters
irf_inf_ye = temp(1:hor+1);         % year end inflaion response in ppts

%% IRFs of interest

irf_gdp_ext  = IRF_ER(3,:);
irf_inf_ext  = irf_inf_ye;
irf_rate_ext = IRF_ER(5,:);
irf_rtwi_ext = IRF_ER(6,:);


%% Extended-sample point-estimate diagnostics

[gdp_trough_ext,gdp_ind_ext] = min(irf_gdp_ext);
[gdp_peak_ext,gdp_peak_ind_ext] = max(irf_gdp_ext);
gdp_horizon_ext = gdp_ind_ext - 1;


[inf_trough_ext,inf_ind_ext] = min(irf_inf_ext);
inf_horizon_ext = inf_ind_ext - 1;

[rate_trough_ext,rate_ind_ext] = min(irf_rate_ext);
rate_horizon_ext = rate_ind_ext - 1;

fprintf('\nExtended-sample exchange-rate shock diagnostics:\n')

fprintf('RTWI impact response: %.3f percent\n', ...
    irf_rtwi_ext(1))

fprintf('GDP trough: %.3f percent at horizon %d quarters\n', ...
    gdp_trough_ext,gdp_horizon_ext)

fprintf('GDP peak: %.3f percent at horizon %d quarters\n', ...
    gdp_peak_ext,gdp_peak_ind_ext-1)

fprintf('Year-ended inflation trough: %.3f ppt at horizon %d quarters\n', ...
    inf_trough_ext,inf_horizon_ext)

fprintf('Cash-rate trough: %.3f ppt at horizon %d quarters\n', ...
    rate_trough_ext,rate_horizon_ext)

%% Extended-sample FEVD

[~,FEVD_40_ext] = ...
    calculate_IRF_FEVD(phi,B,hor_fevd,'FEVD');


%% A.2 analogue - variance explained by exchange-rate shock

h_A2 = [2 10 20 40];
vars_A2 = [3 4 5];

fevd_A2_ext = NaN(length(vars_A2),length(h_A2));

for ii = 1:length(vars_A2)

    fevd_A2_ext(ii,:) = reshape( ...
        FEVD_40_ext(vars_A2(ii),shock_er,h_A2),1,[]);

end

Table_A2_ext = array2table(fevd_A2_ext, ...
    'VariableNames',{'H2','H10','H20','H40'}, ...
    'RowNames',{'GDP','Inflation','CashRate'});

fprintf('\nExtended-sample A.2 analogue:\n')
disp(Table_A2_ext)


%% A.4 analogue - decomposition of RTWI variance

h_A4 = [1 10 20 40];

fevd_foreign_ext = squeeze( ...
    sum(FEVD_40_ext(6,1:2,h_A4),2));

fevd_domestic_ext = squeeze( ...
    sum(FEVD_40_ext(6,3:5,h_A4),2));

fevd_er_ext = squeeze( ...
    FEVD_40_ext(6,6,h_A4));

fevd_A4_ext = [ ...
    fevd_foreign_ext(:)'; ...
    fevd_domestic_ext(:)'; ...
    fevd_er_ext(:)'];

assert(max(abs(sum(fevd_A4_ext,1)-100)) < 1e-8)

Table_A4_ext = array2table(fevd_A4_ext, ...
    'VariableNames',{'H1','H10','H20','H40'}, ...
    'RowNames',{'Foreign', ...
    'DomesticNonExchangeRate', ...
    'ExchangeRate'});

fprintf('\nExtended-sample A.4 analogue:\n')
disp(Table_A4_ext)

%% Load baseline replication results

baseline_file = ...
    fullfile(folder,'results','baseline.mat');

S = load(baseline_file,'results_repl');

baseline = S.results_repl;

clear S

%% Compare point estimates: baseline vs extended sample

fprintf('\nBaseline versus extended sample:\n')
fprintf('                                 Baseline       Extended\n')

fprintf('GDP trough:                    %8.3f       %8.3f\n', ...
    min(baseline.irf_gdp),gdp_trough_ext)

fprintf('YE inflation trough:           %8.3f       %8.3f\n', ...
    min(baseline.irf_inf),inf_trough_ext)

fprintf('Cash-rate trough:              %8.3f       %8.3f\n', ...
    min(baseline.irf_rate),rate_trough_ext)

%% Compare FEVD: extended minus baseline

diff_A2_ext = fevd_A2_ext - baseline.fevd_A2;
diff_A4_ext = fevd_A4_ext - baseline.fevd_A4;

Table_A2_change = array2table(diff_A2_ext, ...
    'VariableNames',{'H2','H10','H20','H40'}, ...
    'RowNames',{'GDP','Inflation','CashRate'});

Table_A4_change = array2table(diff_A4_ext, ...
    'VariableNames',{'H1','H10','H20','H40'}, ...
    'RowNames',{'Foreign', ...
    'DomesticNonExchangeRate', ...
    'ExchangeRate'});

fprintf('\nChange in A.2 analogue: extended - baseline\n')
disp(Table_A2_change)

fprintf('\nChange in A.4 analogue: extended - baseline\n')
disp(Table_A4_change)

%% Compare baseline and extended-sample IRFs

h = 0:hor;

figure;

subplot(2,2,1)
plot(h,baseline.irf_gdp,'--k','LineWidth',1.5)
hold on
plot(h,irf_gdp_ext,'-k','LineWidth',2)
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('GDP')
xlabel('Quarters')
ylabel('%')
legend('1985Q1-2013Q2','1985Q1-2026Q2','Location','best')
grid on

subplot(2,2,2)
plot(h,baseline.irf_inf,'--k','LineWidth',1.5)
hold on
plot(h,irf_inf_ext,'-k','LineWidth',2)
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('Year-ended inflation')
xlabel('Quarters')
ylabel('ppt')
grid on

subplot(2,2,3)
plot(h,baseline.irf_rate,'--k','LineWidth',1.5)
hold on
plot(h,irf_rate_ext,'-k','LineWidth',2)
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('Cash rate')
xlabel('Quarters')
ylabel('ppt')
grid on

subplot(2,2,4)
plot(h,baseline.irf_rtwi,'--k','LineWidth',1.5)
hold on
plot(h,irf_rtwi_ext,'-k','LineWidth',2)
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('Real TWI')
xlabel('Quarters')
ylabel('%')
grid on

%% Save extended-sample point-estimate results

results_ext = struct();

results_ext.sample_start = 1985;
results_ext.sample_end = 2026 + 1/4;

results_ext.p = p;
results_ext.n_foreign = n_foreign;
results_ext.Variables = Variables;

results_ext.phi = phi;
results_ext.gamma = gamma;
results_ext.SIGMA = SIGMA;

results_ext.irf_gdp = irf_gdp_ext;
results_ext.irf_inf = irf_inf_ext;
results_ext.irf_rate = irf_rate_ext;
results_ext.irf_rtwi = irf_rtwi_ext;

results_ext.FEVD_40 = FEVD_40_ext;
results_ext.fevd_A2 = fevd_A2_ext;
results_ext.fevd_A4 = fevd_A4_ext;

results_ext.diff_A2_from_baseline = diff_A2_ext;
results_ext.diff_A4_from_baseline = diff_A4_ext;

save(fullfile(folder,'results','ext_sample_2.full.mat'), ...
    'results_ext')


%% DIAGNOSTICS ==== Inspect large reduced-form residuals

% Standardise each equation's residuals
e_std = e ./ std(e,0,1);

% Largest absolute standardised residual across equations each quarter
max_abs_resid = max(abs(e_std),[],2);

% Show observations from 2019 onward
check_ind = yearlab_ext >= 2019;

resid_check = table( ...
    yearlab_ext(check_ind), ...
    max_abs_resid(check_ind), ...
    'VariableNames',{'Date','MaxAbsStdResidual'});

disp(resid_check)

outlier_ind = max_abs_resid >= 2.5;

outlier_table = table( ...
    yearlab_ext(outlier_ind), ...
    max_abs_resid(outlier_ind), ...
    'VariableNames',{'Date','MaxAbsStdResidual'});

fprintf('\nPotential residual outliers:\n')
disp(outlier_table)

%% Residual outliers by equation

for jj = 1:N

    fprintf('\n%s residuals > 2.5 standard deviations:\n', ...
        Variables{jj})

    ind = abs(e_std(:,jj)) >= 2.5;

    disp(table( ...
        yearlab_ext(ind), ...
        e_std(ind,jj), ...
        'VariableNames',{'Date','StdResidual'}))

end




%% Which equation drives each large residual?

[max_abs_resid, max_eq] = max(abs(e_std),[],2);

check_dates = [2021+2/4; 2021+3/4; 2022];

for ii = 1:length(check_dates)

    ind = yearlab_ext == check_dates(ii);

    fprintf('\nDate: %.2f\n',yearlab_ext(ind))
    fprintf('Largest residual: %s, %.3f std dev\n', ...
        Variables{max_eq(ind)}, ...
        e_std(ind,max_eq(ind)))

    fprintf('All standardized residuals:\n')
    disp(array2table(e_std(ind,:), ...
        'VariableNames',Variables))

end