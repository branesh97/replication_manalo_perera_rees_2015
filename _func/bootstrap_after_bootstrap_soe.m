function [phi_boot,SIGMA_boot,phi_bias_corrected,bias] = ...
    bootstrap_after_bootstrap_soe(phi,gamma,X,e,d,p,n_foreign,nboot1,nboot2)
% BOOTSTRAP_AFTER_BOOTSTRAP_SOE
%
% Kilian (1998) bootstrap-after-bootstrap adapted to the small-open-
% economy VAR used in the Manalo et al. replication.
%
% This follows the structure of bootstrap_after_bootstrap.m used in
% lectures, with modifications to:
%   (1) impose small-open-economy block exogeneity in every bootstrap
%       re-estimation;
%   (2) include deterministic controls in data generation and estimation;
%   (3) store deterministic coefficients separately from phi;
%   (4) calculate bootstrap covariance matrices using residuals from the
%       re-estimated restricted bootstrap VAR.
%
%
% INPUTS
% phi           (1+p*N) x N VAR coefficient matrix
%               [constant; lag coefficients]
%
% gamma         q x N deterministic-control coefficient matrix
%
% X             T x (1+p*N) VAR regressor matrix from olsvar_soe.m
%
% e             T x N reduced-form residual matrix from olsvar_soe.m
%
% d             T x q matrix of deterministic controls
%
% p             VAR lag order
%
% n_foreign     number of variables in foreign block
%
% nboot1        number of first-stage bootstrap replications
%               default = 1000
%
% nboot2        number of second-stage bootstrap replications
%               default = 1000
%
%
% OUTPUTS
% phi_boot              (1+p*N) x N x nboot2
%                       second-stage bootstrap VAR coefficients
%
% SIGMA_boot            N x N x nboot2
%                       second-stage bootstrap covariance matrices
%
% phi_bias_corrected    bias-corrected VAR coefficients used to
%                       generate second-stage bootstrap samples
%
% bias                  estimated bias adjustment, defined as in the
%                       lecture code:
%                       phi - mean(first-stage bootstrap phi)


%% Preliminaries

[T,N] = size(e);
K = size(phi,1);
q = size(d,2);

if nargin < 8 || isempty(nboot1)
    nboot1 = 1000;
end

if nargin < 9 || isempty(nboot2)
    nboot2 = 1000;
end

% Basic checks
if K ~= 1 + p*N
    error('Size of phi is inconsistent with N and p.')
end

if size(X,1) ~= T || size(X,2) ~= K
    error('Dimensions of X are inconsistent with phi/e.')
end

if size(d,1) ~= T
    error('d must have the same number of rows as X and e.')
end

if size(gamma,1) ~= q || size(gamma,2) ~= N
    error('Dimensions of gamma are inconsistent with d and phi.')
end

if n_foreign >= N
    error('n_foreign must be smaller than total number of variables.')
end


%% Columns of X permitted in foreign equations
%
% X is ordered:
% constant,
% lag 1: variables 1,...,N,
% lag 2: variables 1,...,N,
% ...
%
% Foreign equations retain only lags of variables 1:n_foreign.

foreign_ind = 1;

for ll = 1:p

    first_col = 2 + (ll-1)*N;

    foreign_ind = [foreign_ind ...
                   first_col:first_col+n_foreign-1];

end


%% ===============================================================
% FIRST-STAGE BOOTSTRAP
% Estimate finite-sample bias in phi
% ===============================================================

phi_stage_1 = zeros(K,N);

for jj = 1:nboot1

    %% Generate bootstrap sample

    [X_boot,Y_boot] = generate_bootstrap_sample( ...
        phi,gamma,X,e,d,N,T,K);


    %% Re-estimate bootstrap sample imposing SOE restrictions

    [phi_star,~,~,~] = estimate_bootstrap_soe( ...
        X_boot,Y_boot,d,foreign_ind,K,N,n_foreign);


    %% Accumulate first-stage phi estimates

    phi_stage_1 = phi_stage_1 + phi_star;

end


%% Bias correction
%
% Same convention as lecture bootstrap_after_bootstrap.m:
%
% bias = phi - E*(phi_star)
% phi_bias_corrected = phi + bias

bias = phi - phi_stage_1/nboot1;

phi_bias_corrected = phi + bias;


%% Check stability of bias-corrected VAR

Companion = [phi_bias_corrected(2:end,:)'; ...
             eye(N*(p-1)) zeros(N*(p-1),N)];

EIG = max(abs(eig(Companion)));


% If full bias correction produces a non-stationary VAR,
% shrink the correction towards the original estimate.
%
% This follows the logic of the lecture bootstrap code.

delta = 1;

while EIG >= 1

    delta = 0.99*delta;

    phi_bias_corrected = phi + delta*bias;

    Companion = [phi_bias_corrected(2:end,:)'; ...
                 eye(N*(p-1)) zeros(N*(p-1),N)];

    EIG = max(abs(eig(Companion)));

    if delta < 1e-8
        warning(['Could not obtain a stable bias-corrected VAR. ', ...
                 'Using original phi.'])

        phi_bias_corrected = phi;
        break
    end

end


if delta < 1
    fprintf('Bias correction stability adjustment: delta = %.4f\n',delta)
end


%% ===============================================================
% SECOND-STAGE BOOTSTRAP
% Construct sampling distribution around bias-corrected VAR
% ===============================================================

phi_boot = NaN(K,N,nboot2);
SIGMA_boot = NaN(N,N,nboot2);

for jj = 1:nboot2

    %% Generate bootstrap sample using bias-corrected VAR dynamics
    %
    % gamma remains at its original point estimate when generating
    % bootstrap samples. It is nevertheless re-estimated below in
    % every bootstrap replication.

    [X_boot,Y_boot] = generate_bootstrap_sample( ...
        phi_bias_corrected,gamma,X,e,d,N,T,K);


    %% Re-estimate restricted bootstrap VAR

    [phi_star,~,e_star,SIGMA_star] = estimate_bootstrap_soe( ...
        X_boot,Y_boot,d,foreign_ind,K,N,n_foreign);


    %% Store bootstrap estimates

    phi_boot(:,:,jj) = phi_star;

    SIGMA_boot(:,:,jj) = SIGMA_star;

end


end



%% ===============================================================
% LOCAL FUNCTION:
% Generate one artificial bootstrap sample
% ===============================================================

function [X_boot,Y_boot] = generate_bootstrap_sample( ...
    phi_gen,gamma_gen,X,e,d,N,T,K)

% This follows the recursive structure of the lecturer's
% bootstrap_after_bootstrap.m.

X_boot = NaN(T,K);
Y_boot = NaN(T,N);


% Initial conditions:
% same convention as lecture bootstrap code
X_boot(1,:) = X(T,:);


for ii = 1:T

    % Draw an entire reduced-form residual vector.
    % Drawing by row preserves contemporaneous covariance across equations.
    draw = randi(T);


    % Artificial observation
    Y_boot(ii,:) = ...
        X_boot(ii,:)*phi_gen ...
        + d(ii,:)*gamma_gen ...
        + e(draw,:);


    % Construct next observation's lag vector
    if ii < T

        X_boot(ii+1,:) = ...
            [1 Y_boot(ii,:) X_boot(ii,2:end-N)];

    end

end


end



%% ===============================================================
% LOCAL FUNCTION:
% Re-estimate one bootstrap sample under SOE restrictions
% ===============================================================

function [phi_star,gamma_star,e_star,SIGMA_star] = ...
    estimate_bootstrap_soe( ...
    X_boot,Y_boot,d,foreign_ind,K,N,n_foreign)


T = size(Y_boot,1);
q = size(d,2);


%% Regressor matrices

% Foreign equations:
% constant + lags of foreign variables only + deterministic controls
X_foreign = [X_boot(:,foreign_ind) d];

% Domestic equations:
% constant + lags of all variables + deterministic controls
X_domestic = [X_boot d];


%% Storage

phi_star = zeros(K,N);
gamma_star = zeros(q,N);
e_star = NaN(T,N);


%% Foreign equations

for eq = 1:n_foreign

    beta = X_foreign \ Y_boot(:,eq);

    % Constant + permitted foreign lag coefficients
    phi_star(foreign_ind,eq) = ...
        beta(1:length(foreign_ind));

    % Deterministic coefficients
    if q > 0

        gamma_star(:,eq) = ...
            beta(length(foreign_ind)+1:end);

    end

    % Restricted-equation residual
    e_star(:,eq) = ...
        Y_boot(:,eq) - X_foreign*beta;

end


%% Domestic equations

for eq = n_foreign+1:N

    beta = X_domestic \ Y_boot(:,eq);

    % Constant + all VAR lag coefficients
    phi_star(:,eq) = beta(1:K);

    % Deterministic coefficients
    if q > 0

        gamma_star(:,eq) = beta(K+1:end);

    end

    % Residual
    e_star(:,eq) = ...
        Y_boot(:,eq) - X_domestic*beta;

end


%% Reduced-form covariance matrix
%
% Use a common degrees-of-freedom adjustment based on the full domestic
% equation, as in our baseline olsvar_soe.m.

K_full = size(X_domestic,2);

SIGMA_star = (e_star'*e_star)/(T-K_full);


end