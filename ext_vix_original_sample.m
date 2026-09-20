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

%% VIX matched sample: 1990Q1 to 2013Q2

sample_vix = (yearlab >= 1990) & (yearlab <= 2013 + 1/4);

raw_vix = data_raw(sample_vix,:);
yearlab_vix = yearlab(sample_vix);

% All seven variables must be available
assert(~any(isnan(raw_vix(:))))

T = size(raw_vix,1);

% 1990Q1 to 2013Q2
assert(T == 94) % check

%% Model settings
% Extended sample with pandemic controls extension +
% vix for global risk

% therefore two models (WITHOUT VIX and WITH VIX) on same sample

hor = 24;               % IRF horizon: 24 quarters
hor_fevd = 40;

nboot1 = 1000;          % for double bootstrap
nboot2 = 1000;          % for double bootstrap

rng_seed = 12345;       % Fixed seed, replication files produce same bands


assert(T == 94) % check

%% Transform variables

vix = raw_vix(:,7);            % level, quarterly average

log_us_gdp  = 100*log(raw_vix(:,1));
log_tot     = 100*log(raw_vix(:,2));
log_aus_gdp = 100*log(raw_vix(:,3));
log_rtwi    = 100*log(raw_vix(:,6));

trim_inf  = raw_vix(:,4);
cash_rate = raw_vix(:,5);

%% Quadratic detrending of GDP
% IMPORTANT: trends re-estimated over the full extended sample w VIX data

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

% Deterministic controls
% Original controls + pandemic dummies

% (1) Inflation-targeting regime dummy
d_it = double(yearlab_vix >= 1993);

% (2) Individual GFC quarter pulse dummies
d_2008q4 = double(yearlab_vix == 2008 + 3/4);
d_2009q1 = double(yearlab_vix == 2009);
d_2009q2 = double(yearlab_vix == 2009 + 1/4);
d_2009q3 = double(yearlab_vix == 2009 + 2/4);

d = [d_it ...
     d_2008q4 d_2009q1 d_2009q2 d_2009q3]; %...
     %d_2020q2 d_2021q3]; % for neg shock only

%% Matched-sample model WITHOUT VIX

y_noVIX = [us_gdp_gap ...
    log_tot ...
    aus_gdp_gap ...
    trim_inf ...
    cash_rate ...
    log_rtwi];

N_noVIX = size(y_noVIX,2);

p = 2;
n_foreign_noVIX = 2;
shock_er_noVIX = 6;

[phi_noVIX,gamma_noVIX,SIGMA_noVIX,X_noVIX,e_noVIX] = ...
    olsvar_soe(y_noVIX,p,d,n_foreign_noVIX);

B_noVIX = chol(SIGMA_noVIX,'lower');

B_N_noVIX = B_noVIX;
B_N_noVIX(:,shock_er_noVIX) = ...
    10*B_noVIX(:,shock_er_noVIX) / ...
    B_noVIX(shock_er_noVIX,shock_er_noVIX);

IRF_N_noVIX = calculate_IRF_FEVD(phi_noVIX,B_N_noVIX,hor);

ind_er_noVIX = ...
    (shock_er_noVIX-1)*N_noVIX + (1:N_noVIX);

IRF_ER_noVIX = IRF_N_noVIX(ind_er_noVIX,:);

irf_gdp_noVIX = IRF_ER_noVIX(3,:);

temp = conv(IRF_ER_noVIX(4,:),ones(1,4));
irf_inf_noVIX = temp(1:hor+1);

irf_rate_noVIX = IRF_ER_noVIX(5,:);
irf_rtwi_noVIX = IRF_ER_noVIX(6,:);

%% VIX-augmented model

y_vix = [vix ...
    us_gdp_gap ...
    log_tot ...
    aus_gdp_gap ...
    trim_inf ...
    cash_rate ...
    log_rtwi];

Variables_vix = {'VIX', ...
    'US Real GDP', ...
    'Terms of Trade', ...
    'Australian Real GDP', ...
    'Trimmed Mean Inflation', ...
    'Cash Rate', ...
    'Real TWI'};

N_vix = size(y_vix,2);

n_foreign_vix = 3;
shock_er_vix = 7;

[phi_vix,gamma_vix,SIGMA_vix,X_vix,e_vix] = ...
    olsvar_soe(y_vix,p,d,n_foreign_vix);

B_vix = chol(SIGMA_vix,'lower');

%% Matched-sample no-VIX FEVD

[~,FEVD_40_noVIX] = ...
    calculate_IRF_FEVD(phi_noVIX,B_noVIX,hor_fevd,'FEVD');

h_A4 = [1 10 20 40];

% RTWI = variable 6 in no-VIX model

fevd_foreign_noVIX = squeeze( ...
    sum(FEVD_40_noVIX(6,1:2,h_A4),2));

fevd_domestic_noVIX = squeeze( ...
    sum(FEVD_40_noVIX(6,3:5,h_A4),2));

fevd_er_noVIX = squeeze( ...
    FEVD_40_noVIX(6,6,h_A4));

fevd_A4_noVIX = [ ...
    fevd_foreign_noVIX(:)'; ...
    fevd_domestic_noVIX(:)'; ...
    fevd_er_noVIX(:)'];

Table_A4_noVIX = array2table(fevd_A4_noVIX, ...
    'VariableNames',{'H1','H10','H20','H40'}, ...
    'RowNames',{'Foreign', ...
    'DomesticNonExchangeRate', ...
    'ExchangeRate'});

fprintf('\nRTWI FEVD - matched sample without VIX:\n')
disp(Table_A4_noVIX)

assert(max(abs(sum(fevd_A4_noVIX,1)-100)) < 1e-8)

B_N_vix = B_vix;
B_N_vix(:,shock_er_vix) = ...
    10*B_vix(:,shock_er_vix) / ...
    B_vix(shock_er_vix,shock_er_vix);

IRF_N_vix = calculate_IRF_FEVD(phi_vix,B_N_vix,hor);

ind_er_vix = ...
    (shock_er_vix-1)*N_vix + (1:N_vix);

IRF_ER_vix = IRF_N_vix(ind_er_vix,:);

irf_gdp_vix = IRF_ER_vix(4,:);

temp = conv(IRF_ER_vix(5,:),ones(1,4));
irf_inf_vix = temp(1:hor+1);

irf_rate_vix = IRF_ER_vix(6,:);
irf_rtwi_vix = IRF_ER_vix(7,:);

%% Compare exchange-rate IRFs: same sample, without vs with VIX

h = 0:hor;

figure;

subplot(2,2,1)
plot(h,irf_gdp_noVIX,'--k','LineWidth',1.5)
hold on
plot(h,irf_gdp_vix,'-k','LineWidth',2)
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('GDP')
xlabel('Quarters')
ylabel('%')
legend('No VIX','With VIX','Location','best')
grid on

subplot(2,2,2)
plot(h,irf_inf_noVIX,'--k','LineWidth',1.5)
hold on
plot(h,irf_inf_vix,'-k','LineWidth',2)
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('Year-ended inflation')
xlabel('Quarters')
ylabel('ppt')
grid on

subplot(2,2,3)
plot(h,irf_rate_noVIX,'--k','LineWidth',1.5)
hold on
plot(h,irf_rate_vix,'-k','LineWidth',2)
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('Cash rate')
xlabel('Quarters')
ylabel('ppt')
grid on

subplot(2,2,4)
plot(h,irf_rtwi_noVIX,'--k','LineWidth',1.5)
hold on
plot(h,irf_rtwi_vix,'-k','LineWidth',2)
plot([0 hor],[0 0],'-k','LineWidth',0.75)
title('Real TWI')
xlabel('Quarters')
ylabel('%')
grid on

%% VIX-model FEVD

[~,FEVD_40_vix] = ...
    calculate_IRF_FEVD(phi_vix,B_vix,hor_fevd,'FEVD');

h_A4 = [1 10 20 40];

% RTWI is variable 7

fevd_vixshock = squeeze( ...
    FEVD_40_vix(7,1,h_A4));

fevd_foreignmacro = squeeze( ...
    sum(FEVD_40_vix(7,2:3,h_A4),2));

fevd_domestic_vix = squeeze( ...
    sum(FEVD_40_vix(7,4:6,h_A4),2));

fevd_er_vix = squeeze( ...
    FEVD_40_vix(7,7,h_A4));

fevd_A4_vix = [ ...
    fevd_vixshock(:)'; ...
    fevd_foreignmacro(:)'; ...
    fevd_domestic_vix(:)'; ...
    fevd_er_vix(:)'];

Table_A4_vix = array2table(fevd_A4_vix, ...
    'VariableNames',{'H1','H10','H20','H40'}, ...
    'RowNames',{'VIX', ...
    'OtherForeign', ...
    'DomesticNonExchangeRate', ...
    'ExchangeRate'});

fprintf('\nRTWI FEVD with VIX:\n')
disp(Table_A4_vix)

assert(max(abs(sum(fevd_A4_vix,1)-100)) < 1e-8)

foreign_total_vix = ...
    fevd_vixshock + fevd_foreignmacro;

fevd_A4_vix_grouped = [ ...
    foreign_total_vix(:)'; ...
    fevd_domestic_vix(:)'; ...
    fevd_er_vix(:)'];

%% Effect of adding VIX on matched-sample RTWI FEVD

diff_A4_vix = fevd_A4_vix_grouped - fevd_A4_noVIX;

Table_A4_vix_change = array2table(diff_A4_vix, ...
    'VariableNames',{'H1','H10','H20','H40'}, ...
    'RowNames',{'Foreign', ...
    'DomesticNonExchangeRate', ...
    'ExchangeRate'});

fprintf('\nChange from adding VIX: with VIX - without VIX\n')
disp(Table_A4_vix_change)