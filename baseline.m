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

%% Original replication sample: 1985Q1 to 2013Q2

sample_repl = (yearlab >= 1985) & (yearlab <= 2013 + 1/4);

raw_repl = data_raw(sample_repl,1:6);   % VIX not used in baseline model
assert(~any(isnan(raw_repl(:))))        % check
yearlab_repl = yearlab(sample_repl);

T = size(raw_repl,1);
assert(T == 114) % check

%% Model settings

p = 2;                  % Manalo et al. use two lags
n_foreign = 2;          % US GDP and terms of trade
shock_er = 6;           % Real TWI is ordered last

hor = 24;               % IRF horizon: 24 quarters
hor_fevd = 40;

nboot1 = 1000;          % for double bootstrap
nboot2 = 1000;          % for double bootstrap

rng_seed = 12345;       % Fixed seed, replication files produce same bands


%% Transform variables

% Log variables
log_us_gdp  = 100*log(raw_repl(:,1));
log_tot     = 100*log(raw_repl(:,2));
log_aus_gdp = 100*log(raw_repl(:,3));
log_rtwi    = 100*log(raw_repl(:,6));

% Already in required units
trim_inf  = raw_repl(:,4);      % quarterly percentage change
cash_rate = raw_repl(:,5);      % percent, level

%% Quadratic detrending of GDP

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

% Deterministic controls used by Manalo et al. (2015)

% (1) Inflation-targeting regime dummy
d_it = double(yearlab_repl >= 1993);

% (2) Individual GFC quarter pulse dummies
d_2008q4 = double(yearlab_repl == 2008 + 3/4);
d_2009q1 = double(yearlab_repl == 2009);
d_2009q2 = double(yearlab_repl == 2009 + 1/4);
d_2009q3 = double(yearlab_repl == 2009 + 2/4);

d = [d_it d_2008q4 d_2009q1 d_2009q2 d_2009q3];

%% Estimate aggregate VAR

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

%% Figure 1 variables

irf_gdp  = IRF_ER(3,:);       % percent deviation from quadratic trend
irf_inf  = irf_inf_ye;        % year-ended inflation, percentage points
irf_rate = IRF_ER(5,:);       % percentage points
irf_rtwi = IRF_ER(6,:);       % percent

%% Point-estimate diagnostics

% GDP trough
[gdp_trough,gdp_ind] = min(irf_gdp);
gdp_horizon = gdp_ind - 1;

% Year-ended inflation trough
[inf_trough,inf_ind] = min(irf_inf);
inf_horizon = inf_ind - 1;

% Cash-rate trough
[rate_trough,rate_ind] = min(irf_rate);
rate_horizon = rate_ind - 1;

fprintf('\nExchange-rate shock diagnostics:\n')
fprintf('RTWI impact response: %.3f percent\n',irf_rtwi(1))
fprintf('GDP trough: %.3f percent at horizon %d quarters\n', ...
    gdp_trough,gdp_horizon)
fprintf('Year-ended inflation trough: %.3f ppt at horizon %d quarters\n', ...
    inf_trough,inf_horizon)
fprintf('Cash-rate trough: %.3f ppt at horizon %d quarters\n', ...
    rate_trough,rate_horizon)

%% Kilian (1998) bootstrap-after-bootstrap

rng(rng_seed)

[phi_boot,SIGMA_boot,phi_bias_corrected,bias] = ...
    bootstrap_after_bootstrap_soe( ...
    phi,gamma,X,e,d,p,n_foreign,nboot1,nboot2);

%% Bootstrap impulse response functions - all draws

IRF_N_boot = NaN(N^2,hor+1,nboot2);

for jj = 1:nboot2

    % Cholesky identification
    B_boot = chol(SIGMA_boot(:,:,jj),'lower');

    % Normalize exchange-rate shock to 10 percent appreciation
    B_N_boot = B_boot;
    B_N_boot(:,shock_er) = ...
        10*B_boot(:,shock_er)./B_boot(shock_er,shock_er);

    % Bootstrap IRFs
    IRF_N_boot(:,:,jj) = ...
        calculate_IRF_FEVD(phi_boot(:,:,jj),B_N_boot,hor);

end

%% Extract bootstrap exchange-rate shock responses

IRF_ER_boot = IRF_N_boot(ind_er,:,:);

irf_gdp_boot   = squeeze(IRF_ER_boot(3,:,:));
irf_inf_q_boot = squeeze(IRF_ER_boot(4,:,:));
irf_rate_boot  = squeeze(IRF_ER_boot(5,:,:));
irf_rtwi_boot  = squeeze(IRF_ER_boot(6,:,:));

%% Quarterly inflation to year-ended inflation

irf_inf_ye_boot = NaN(hor+1,nboot2);

for jj = 1:nboot2

    temp_boot = conv(irf_inf_q_boot(:,jj)',ones(1,4));

    irf_inf_ye_boot(:,jj) = ...
        temp_boot(1:hor+1)';

end

%% Bootstrap confidence intervals - main specification

CI_quantiles = [0.025 0.075 0.925 0.975];

gdp_CI  = quantile(irf_gdp_boot,CI_quantiles,2);
inf_CI  = quantile(irf_inf_ye_boot,CI_quantiles,2);
rate_CI = quantile(irf_rate_boot,CI_quantiles,2);
rtwi_CI = quantile(irf_rtwi_boot,CI_quantiles,2);

%% Figure 1 replication with bootstrap confidence intervals

h = 0:hor;

figure;

% GDP
subplot(2,2,1)
plot(h,irf_gdp,'-k','LineWidth',2)
hold on
plot(h,gdp_CI(:,2:3),'--b','LineWidth',1.2)      % 85%
plot(h,gdp_CI(:,[1 4]),'--r','LineWidth',1.2)    % 95%
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('GDP')
xlabel('Quarters')
ylabel('%')
xlim([0 hor])
grid on

% Year-ended inflation
subplot(2,2,2)
plot(h,irf_inf,'-k','LineWidth',2)
hold on
plot(h,inf_CI(:,2:3),'--b','LineWidth',1.2)
plot(h,inf_CI(:,[1 4]),'--r','LineWidth',1.2)
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('Year-ended inflation')
xlabel('Quarters')
ylabel('ppt')
xlim([0 hor])
grid on

% Cash rate
subplot(2,2,3)
plot(h,irf_rate,'-k','LineWidth',2)
hold on
plot(h,rate_CI(:,2:3),'--b','LineWidth',1.2)
plot(h,rate_CI(:,[1 4]),'--r','LineWidth',1.2)
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('Cash rate')
xlabel('Quarters')
ylabel('ppt')
xlim([0 hor])
grid on

% Real TWI
subplot(2,2,4)
plot(h,irf_rtwi,'-k','LineWidth',2)
hold on
plot(h,rtwi_CI(:,2:3),'--b','LineWidth',1.2)
plot(h,rtwi_CI(:,[1 4]),'--r','LineWidth',1.2)
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('Real TWI')
xlabel('Quarters')
ylabel('%')
xlim([0 hor])
grid on

%% FEVD - aggregate model
% Replicate Tables A.2 and A.4 of Manalo et al. (2015)

% IMPORTANT: use the original Cholesky impact matrix B,
% NOT the 10-percent-normalized B_N.
%
% FEVD assumes the structural shocks have unit variance.
[~,FEVD_40] = calculate_IRF_FEVD(phi,B,hor_fevd,'FEVD');

%% Table A.2
% Percentage of variance explained by exchange-rate shocks

h_A2 = [2 10 20 40];
vars_A2 = [3 4 5];

fevd_A2 = NaN(length(vars_A2),length(h_A2));

for ii = 1:length(vars_A2)

    fevd_A2(ii,:) = reshape( ...
        FEVD_40(vars_A2(ii),shock_er,h_A2),1,[]);

end

Table_A2 = array2table(fevd_A2, ...
    'VariableNames',{'H2','H10','H20','H40'}, ...
    'RowNames',{'GDP','Inflation','CashRate'});

fprintf('\nTable A.2 replication:\n')
disp(Table_A2)

%% Compare Table A.2 with published values

paper_A2 = [ ...
    0.2  2.6  2.1  1.5;   % GDP
    1.0  3.6  3.9  3.3;   % Inflation
    0.1  1.9  3.8  3.2];  % Cash rate

diff_A2 = fevd_A2 - paper_A2;

Table_A2_diff = array2table(diff_A2, ...
    'VariableNames',{'H2','H10','H20','H40'}, ...
    'RowNames',{'GDP','Inflation','CashRate'});

fprintf('\nDifference from published Table A.2 (replication - paper):\n')
disp(Table_A2_diff)

%% Table A.4
% Decomposition of real exchange rate variance

h_A4 = [1 10 20 40];

% Variable 6 = real TWI

fevd_foreign = squeeze( ...
    sum(FEVD_40(6,1:2,h_A4),2));

fevd_domestic = squeeze( ...
    sum(FEVD_40(6,3:5,h_A4),2));

fevd_er = squeeze( ...
    FEVD_40(6,6,h_A4));

% Put into rows
fevd_A4 = [ ...
    fevd_foreign(:)'; ...
    fevd_domestic(:)'; ...
    fevd_er(:)'];
assert(max(abs(sum(fevd_A4,1)-100)) < 1e-8) % check

Table_A4 = array2table(fevd_A4, ...
    'VariableNames',{'H1','H10','H20','H40'}, ...
    'RowNames',{'Foreign', ...
    'DomesticNonExchangeRate', ...
    'ExchangeRate'});

fprintf('\nTable A.4 replication:\n')
disp(Table_A4)

%% Compare Table A.4 with published values

paper_A4 = [ ...
    11.3  60.1  82.0  92.0;   % Foreign
    2.6  11.8   5.4   2.4;   % Domestic non-exchange-rate
    86.1  28.1  12.6   5.6];  % Exchange rate

diff_A4 = fevd_A4 - paper_A4;

Table_A4_diff = array2table(diff_A4, ...
    'VariableNames',{'H1','H10','H20','H40'}, ...
    'RowNames',{'Foreign', ...
    'DomesticNonExchangeRate', ...
    'ExchangeRate'});

fprintf('\nDifference from published Table A.4 (replication - paper):\n')
disp(Table_A4_diff)

%% Save baseline replication results

results_repl = struct();

% Model/sample information
results_repl.sample_start = 1985;
results_repl.sample_end = 2013 + 1/4;
results_repl.p = p;
results_repl.n_foreign = n_foreign;
results_repl.shock_er = shock_er;
results_repl.rng_seed = rng_seed;
results_repl.Variables = Variables;

% VAR estimates
results_repl.phi = phi;
results_repl.gamma = gamma;
results_repl.SIGMA = SIGMA;

% Figure 1
results_repl.irf_gdp  = irf_gdp;
results_repl.irf_inf  = irf_inf;
results_repl.irf_rate = irf_rate;
results_repl.irf_rtwi = irf_rtwi;

results_repl.gdp_CI  = gdp_CI;
results_repl.inf_CI  = inf_CI;
results_repl.rate_CI = rate_CI;
results_repl.rtwi_CI = rtwi_CI;

% FEVD
results_repl.FEVD_40 = FEVD_40;

results_repl.fevd_A2 = fevd_A2;
results_repl.paper_A2 = paper_A2;
results_repl.diff_A2 = diff_A2;

results_repl.fevd_A4 = fevd_A4;
results_repl.paper_A4 = paper_A4;
results_repl.diff_A4 = diff_A4;

% Bootstrap information
results_repl.phi_bias_corrected = phi_bias_corrected;
results_repl.bias = bias;

save(fullfile(folder,'results','baseline.mat'), 'results_repl')