function [phi,gamma,SIGMA,X,e] = olsvar_soe(y,p,D,n_foreign)
% Based on: 
% OLSVAR by
% Benjamin Wong
% Monash University
% Edited July 2019
% Estimates a VAR with a constant using least squares with options of
% various methods of bootstrapping
%
% OLSVAR_SOE
% Estimates a VAR with small-open-economy block exogeneity.
% Branesh Prakash (created September 2026)
%
% The first n_foreign variables are treated as the foreign block:
%   - foreign equations contain only lags of foreign variables
%   - domestic equations contain lags of all variables
%
% Deterministic/exogenous controls D are included in all equations,
% but their coefficients are stored separately in gamma so that phi
% retains the format required by calculate_IRF_FEVD.m.
%
% Mean backcasting of the first p observations follows olsvar.m.
%
% INPUTS
% y             T x N matrix of endogenous variables
% p             VAR lag order
% D             T x q matrix of deterministic controls
%               (can be [] if none)
% n_foreign     number of variables in foreign block
%
% OUTPUTS
% phi           (1+p*N) x N matrix:
%               constant + VAR lag coefficients
%               restricted coefficients are stored as zeros
% gamma         q x N matrix of coefficients on D
% SIGMA         N x N covariance matrix of reduced-form residuals
% X             T x (1+p*N) unrestricted VAR regressor matrix
% e             T x N matrix of reduced-form residuals

%% Preliminaries

[T,N] = size(y);

if nargin < 3 || isempty(D)
    D = zeros(T,0);
end

if nargin < 4
    n_foreign = 2;
end

% Checks
if size(D,1) ~= T
    error('D must have the same number of rows as y.')
end

if n_foreign >= N
    error('n_foreign must be smaller than the total number of variables.')
end

q = size(D,2);

%% Backcast first p observations using sample mean
% Same convention as olsvar.m

y_bc = [repmat(mean(y,1),p,1); y];


%% Construct dependent-variable and lag matrices

Y = y_bc(p+1:end,:);

% Full VAR regressor matrix:
% constant, lag 1 of all N variables, ..., lag p of all N variables
X = ones(T,1);

for ii = 1:p
    Z = y_bc(p+1-ii:end-ii,:);
    X = [X Z];
end

%% Identify columns of X allowed in foreign equations

% First column is the constant
foreign_ind = 1;

% Within each lag block, retain only the first n_foreign variables
for ii = 1:p
    first_col = 2 + (ii-1)*N;
    foreign_ind = [foreign_ind ...
        first_col:first_col+n_foreign-1];
end

% Regressor matrices including deterministic controls
X_foreign = [X(:,foreign_ind) D];
X_domestic = [X D];

%% Estimate equations

% Keep phi in standard VAR format expected by calculate_IRF_FEVD
phi = zeros(1+p*N,N);

% Store deterministic coefficients separately
gamma = zeros(q,N);

% Reduced-form residuals
e = NaN(T,N);


% Foreign equations
for jj = 1:n_foreign

    beta = X_foreign \ Y(:,jj);

    % Constant and permitted foreign lag coefficients
    phi(foreign_ind,jj) = beta(1:length(foreign_ind));

    % Deterministic controls
    if q > 0
        gamma(:,jj) = beta(length(foreign_ind)+1:end);
    end

    % Residuals
    e(:,jj) = Y(:,jj) - X_foreign*beta;

end


% Domestic equations
for jj = n_foreign+1:N

    beta = X_domestic \ Y(:,jj);

    % Constant and all VAR lag coefficients
    phi(:,jj) = beta(1:size(X,2));

    % Deterministic controls
    if q > 0
        gamma(:,jj) = beta(size(X,2)+1:end);
    end

    % Residuals
    e(:,jj) = Y(:,jj) - X_domestic*beta;

end


%% Reduced-form residual covariance matrix

% Use a common degrees-of-freedom adjustment so SIGMA remains a valid
% covariance matrix. K_full is the number of regressors in a domestic
% equation, including deterministic controls.
K_full = size(X_domestic,2);

SIGMA = (e'*e)/(T-K_full);


end