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

nboot1 = 1000;          % for double bootstrap
nboot2 = 1000;          % for double bootstrap

rng_seed = 12345;       % Fixed seed, replication files produce same bands


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
% Original controls + pandemic dummies

% (1) Inflation-targeting regime dummy
d_it = double(yearlab_ext >= 1993);

% (2) Individual GFC quarter pulse dummies
d_2008q4 = double(yearlab_ext == 2008 + 3/4);
d_2009q1 = double(yearlab_ext == 2009);
d_2009q2 = double(yearlab_ext == 2009 + 1/4);
d_2009q3 = double(yearlab_ext == 2009 + 2/4);

% (3) Pandemic controls
d_2020q2 = double(yearlab_ext == 2020 + 1/4);
d_2021q3 = double(yearlab_ext == 2021 + 2/4);

d = [d_it ...
     d_2008q4 d_2009q1 d_2009q2 d_2009q3 ...
     d_2020q2 d_2021q3]; % for neg shock only

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

fprintf('\nExtended-sample (with pandemic dummy) exchange-rate shock diagnostics:\n')

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

%% Kilian (1998) bootstrap-after-bootstrap
% Preferred extended sample with pandemic controls

rng(rng_seed)

[phi_boot,SIGMA_boot,phi_bias_corrected,bias] = ...
    bootstrap_after_bootstrap_soe( ...
    phi,gamma,X,e,d,p,n_foreign,nboot1,nboot2);

%% Bootstrap impulse response functions

IRF_N_boot = NaN(N^2,hor+1,nboot2);

for jj = 1:nboot2

    % Cholesky identification
    B_boot = chol(SIGMA_boot(:,:,jj),'lower');

    % Normalize exchange-rate shock to a 10 percent appreciation
    B_N_boot = B_boot;

    B_N_boot(:,shock_er) = ...
        10*B_boot(:,shock_er)./B_boot(shock_er,shock_er);

    % Bootstrap normalized IRFs
    IRF_N_boot(:,:,jj) = ...
        calculate_IRF_FEVD(phi_boot(:,:,jj),B_N_boot,hor);

end

%% Extract bootstrap exchange-rate shock responses

IRF_ER_boot = IRF_N_boot(ind_er,:,:);

irf_gdp_boot   = squeeze(IRF_ER_boot(3,:,:));
irf_inf_q_boot = squeeze(IRF_ER_boot(4,:,:));
irf_rate_boot  = squeeze(IRF_ER_boot(5,:,:));
irf_rtwi_boot  = squeeze(IRF_ER_boot(6,:,:));

%% Convert bootstrap quarterly inflation responses to year-ended inflation

irf_inf_ye_boot = NaN(hor+1,nboot2);

for jj = 1:nboot2

    temp_boot = conv(irf_inf_q_boot(:,jj)',ones(1,4));

    irf_inf_ye_boot(:,jj) = ...
        temp_boot(1:hor+1)';

end

%% Bootstrap confidence intervals

CI_quantiles = [0.025 0.075 0.925 0.975];

gdp_CI_ext  = quantile(irf_gdp_boot,CI_quantiles,2);
inf_CI_ext  = quantile(irf_inf_ye_boot,CI_quantiles,2);
rate_CI_ext = quantile(irf_rate_boot,CI_quantiles,2);
rtwi_CI_ext = quantile(irf_rtwi_boot,CI_quantiles,2);

%% Preferred extended-sample IRFs with bootstrap confidence intervals

h = 0:hor;

figure;

% GDP
subplot(2,2,1)
plot(h,irf_gdp_ext,'-k','LineWidth',2)
hold on
plot(h,gdp_CI_ext(:,2:3),'--b','LineWidth',1.2)      % 85%
plot(h,gdp_CI_ext(:,[1 4]),'--r','LineWidth',1.2)    % 95%
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('GDP')
xlabel('Quarters')
ylabel('%')
xlim([0 hor])
grid on

% Year-ended inflation
subplot(2,2,2)
plot(h,irf_inf_ext,'-k','LineWidth',2)
hold on
plot(h,inf_CI_ext(:,2:3),'--b','LineWidth',1.2)
plot(h,inf_CI_ext(:,[1 4]),'--r','LineWidth',1.2)
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('Year-ended inflation')
xlabel('Quarters')
ylabel('ppt')
xlim([0 hor])
grid on

% Cash rate
subplot(2,2,3)
plot(h,irf_rate_ext,'-k','LineWidth',2)
hold on
plot(h,rate_CI_ext(:,2:3),'--b','LineWidth',1.2)
plot(h,rate_CI_ext(:,[1 4]),'--r','LineWidth',1.2)
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('Cash rate')
xlabel('Quarters')
ylabel('ppt')
xlim([0 hor])
grid on

% Real TWI
subplot(2,2,4)
plot(h,irf_rtwi_ext,'-k','LineWidth',2)
hold on
plot(h,rtwi_CI_ext(:,2:3),'--b','LineWidth',1.2)
plot(h,rtwi_CI_ext(:,[1 4]),'--r','LineWidth',1.2)
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('Real TWI')
xlabel('Quarters')
ylabel('%')
xlim([0 hor])
grid on

%% Bootstrap diagnostics

% Original VAR stability
Companion_original = [phi(2:end,:)'; ...
    eye(N*(p-1)) zeros(N*(p-1),N)];

root_original = max(abs(eig(Companion_original)));

% Bias-corrected VAR stability
Companion_bc = [phi_bias_corrected(2:end,:)'; ...
    eye(N*(p-1)) zeros(N*(p-1),N)];

root_bc = max(abs(eig(Companion_bc)));

fprintf('\nBootstrap diagnostics:\n')
fprintf('Original VAR max root:       %.4f\n',root_original)
fprintf('Bias-corrected VAR max root: %.4f\n',root_bc)
fprintf('Maximum absolute bias:       %.4f\n',max(abs(bias(:))))

fprintf('phi_boot NaNs:   %d\n',sum(isnan(phi_boot(:))))
fprintf('SIGMA_boot NaNs: %d\n',sum(isnan(SIGMA_boot(:))))

%% Stability of second-stage bootstrap draws

unstable_count = 0;
max_root_boot = NaN(nboot2,1);

for jj = 1:nboot2

    Companion_boot = [phi_boot(2:end,:,jj)'; ...
        eye(N*(p-1)) zeros(N*(p-1),N)];

    max_root_boot(jj) = ...
        max(abs(eig(Companion_boot)));

    if max_root_boot(jj) >= 1
        unstable_count = unstable_count + 1;
    end

end

fprintf('\nSecond-stage bootstrap stability:\n')
fprintf('Unstable draws: %d out of %d\n', ...
    unstable_count,nboot2)
fprintf('Median maximum root: %.4f\n', ...
    median(max_root_boot))
fprintf('95th percentile maximum root: %.4f\n', ...
    quantile(max_root_boot,0.95))

fprintf('\nMaximum absolute bootstrap IRFs:\n')
fprintf('GDP:          %.3f\n',max(abs(irf_gdp_boot(:))))
fprintf('YE inflation: %.3f\n',max(abs(irf_inf_ye_boot(:))))
fprintf('Cash rate:    %.3f\n',max(abs(irf_rate_boot(:))))
fprintf('Real TWI:     %.3f\n',max(abs(irf_rtwi_boot(:))))

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

fprintf('\nExtended-sample + pandemic A.2 analogue:\n')
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

fprintf('\nExtended-sample + pandemic A.4 analogue:\n')
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

fprintf('\nChange in A.2 analogue: extended (+ pandemic) - baseline\n')
disp(Table_A2_change)

fprintf('\nChange in A.4 analogue: extended (+ pandemic) - baseline\n')
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

%% Load extended sample (no pandemic) results

extsample_file = ...
    fullfile(folder,'results','ext_sample_2.full.mat');

S = load(extsample_file,'results_ext');

extsample = S.results_ext;

clear S

%% Compare point estimates: extended sample (no pandemic vs pandemic)

fprintf('Extended sample (no pandemic vs pandemic):\n')
fprintf('                                 Ext (no pandemic)       Ext + pandemic\n')

fprintf('GDP trough:                    %8.3f           %8.3f\n', ...
    min(extsample.irf_gdp),gdp_trough_ext)

fprintf('YE inflation trough:           %8.3f           %8.3f\n', ...
    min(extsample.irf_inf),inf_trough_ext)

fprintf('Cash-rate trough:              %8.3f           %8.3f\n', ...
    min(extsample.irf_rate),rate_trough_ext)

%% Compare FEVD: extended minus baseline

diff_A2_ext = fevd_A2_ext - extsample.fevd_A2;
diff_A4_ext = fevd_A4_ext - extsample.fevd_A4;

Table_A2_change = array2table(diff_A2_ext, ...
    'VariableNames',{'H2','H10','H20','H40'}, ...
    'RowNames',{'GDP','Inflation','CashRate'});

Table_A4_change = array2table(diff_A4_ext, ...
    'VariableNames',{'H1','H10','H20','H40'}, ...
    'RowNames',{'Foreign', ...
    'DomesticNonExchangeRate', ...
    'ExchangeRate'});

fprintf('\nChange in A.2 analogue: ext + pandemic - ext (no pandemic)\n')
disp(Table_A2_change)

fprintf('\nChange in A.4 analogue: ext + pandemic - ext (no pandemic)\n')
disp(Table_A4_change)


%% Compare extended-sample (no pandemic vs pandemic) IRFs

h = 0:hor;

figure;

subplot(2,2,1)
plot(h,extsample.irf_gdp,'--k','LineWidth',1.5)
hold on
plot(h,irf_gdp_ext,'-k','LineWidth',2)
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('GDP')
xlabel('Quarters')
ylabel('%')
legend('No pandemic controls','With pandemic dummies','Location','best')
grid on

subplot(2,2,2)
plot(h,extsample.irf_inf,'--k','LineWidth',1.5)
hold on
plot(h,irf_inf_ext,'-k','LineWidth',2)
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('Year-ended inflation')
xlabel('Quarters')
ylabel('ppt')
grid on

subplot(2,2,3)
plot(h,extsample.irf_rate,'--k','LineWidth',1.5)
hold on
plot(h,irf_rate_ext,'-k','LineWidth',2)
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('Cash rate')
xlabel('Quarters')
ylabel('ppt')
grid on

subplot(2,2,4)
plot(h,extsample.irf_rtwi,'--k','LineWidth',1.5)
hold on
plot(h,irf_rtwi_ext,'-k','LineWidth',2)
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('Real TWI')
xlabel('Quarters')
ylabel('%')
grid on

%% Statistical significance horizons

% Columns:
% 1 = lower 95
% 2 = lower 85
% 3 = upper 85
% 4 = upper 95

sig95_gdp  = (gdp_CI_ext(:,1) > 0) | (gdp_CI_ext(:,4) < 0);
sig85_gdp  = (gdp_CI_ext(:,2) > 0) | (gdp_CI_ext(:,3) < 0);

sig95_inf   = (inf_CI_ext(:,1) > 0) | (inf_CI_ext(:,4) < 0);
sig85_inf   = (inf_CI_ext(:,2) > 0) | (inf_CI_ext(:,3) < 0);

sig95_rate  = (rate_CI_ext(:,1) > 0) | (rate_CI_ext(:,4) < 0);
sig85_rate  = (rate_CI_ext(:,2) > 0) | (rate_CI_ext(:,3) < 0);

sig95_rtwi  = (rtwi_CI_ext(:,1) > 0) | (rtwi_CI_ext(:,4) < 0);
sig85_rtwi  = (rtwi_CI_ext(:,2) > 0) | (rtwi_CI_ext(:,3) < 0);

fprintf('\nHorizons significant at 95%%:\n')
fprintf('GDP:       %s\n',mat2str(find(sig95_gdp)'-1))
fprintf('Inflation: %s\n',mat2str(find(sig95_inf)'-1))
fprintf('Cash rate: %s\n',mat2str(find(sig95_rate)'-1))
fprintf('Real TWI:  %s\n',mat2str(find(sig95_rtwi)'-1))

fprintf('\nHorizons significant at 85%%:\n')
fprintf('GDP:       %s\n',mat2str(find(sig85_gdp)'-1))
fprintf('Inflation: %s\n',mat2str(find(sig85_inf)'-1))
fprintf('Cash rate: %s\n',mat2str(find(sig85_rate)'-1))
fprintf('Real TWI:  %s\n',mat2str(find(sig85_rtwi)'-1))

%% Save preferred extended-sample results

results_ext = struct();

% Model/sample information
results_ext.sample_start = 1985;
results_ext.sample_end = 2026 + 1/4;

results_ext.p = p;
results_ext.n_foreign = n_foreign;
results_ext.shock_er = shock_er;
results_ext.rng_seed = rng_seed;
results_ext.Variables = Variables;

% Pandemic controls
results_ext.pandemic_controls = {'2020Q2','2021Q3'};

% VAR estimates
results_ext.phi = phi;
results_ext.gamma = gamma;
results_ext.SIGMA = SIGMA;

% Point IRFs
results_ext.irf_gdp = irf_gdp_ext;
results_ext.irf_inf = irf_inf_ext;
results_ext.irf_rate = irf_rate_ext;
results_ext.irf_rtwi = irf_rtwi_ext;

% Bootstrap confidence intervals
results_ext.gdp_CI = gdp_CI_ext;
results_ext.inf_CI = inf_CI_ext;
results_ext.rate_CI = rate_CI_ext;
results_ext.rtwi_CI = rtwi_CI_ext;

% Bootstrap information
results_ext.phi_bias_corrected = phi_bias_corrected;
results_ext.bias = bias;
results_ext.max_root_boot = max_root_boot;
results_ext.unstable_boot_draws = unstable_count;

% FEVD
results_ext.FEVD_40 = FEVD_40_ext;
results_ext.fevd_A2 = fevd_A2_ext;
results_ext.fevd_A4 = fevd_A4_ext;

save(fullfile(folder,'results', ...
    'ext_sample_3.pandemic_controls.mat'), ...
    'results_ext')
